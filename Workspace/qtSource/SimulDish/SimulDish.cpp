/*Author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it

SimulDish evaluates the flux distribution on the focal plane of a dish on the basis of
the 3D-shape of the dish and the Sun Conic Reflectance of its reflecting surface


   Copyright (C) 2022  Marco Montecchi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include "SimulDish.h"
#include <QtGui>
#include <QMessageBox>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <unistd.h>
#include <qwt_color_map.h>
#include <qwt_plot_spectrogram.h>
#include <qwt_plot_layout.h>
#include <qwt_matrix_raster_data.h>
#include <qwt_scale_widget.h>
#include <qwt_plot_magnifier.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_renderer.h>
#include <qwt_plot_grid.h>
#include <qwt_plot_canvas.h>
#include <qfiledialog.h>
#include <qimagewriter.h>
#include <qwt_plot.h>
#include <qwt_plot_marker.h>
#include <qwt_symbol.h>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia 
Roma - Italy*/

// functions
  double Dr(double r);
  double reflection(double r, double alpha);
  void TrasformXYZ(double yt , double bet);
  double unitVectorReflectance();
  void setBGR(double flux);
 
// known everywhere
  double pig=3.14159265359;
  int S1P2;//Sphere = 1 ; Paraboloid = 2; Experimental(kVISish) = 3
  int daisy;//0 central dish; 1 four identical mirrors arranged like a daisy
  double ds=10.0;  // mm discretization step
  double r1,r2,alpha1,alpha2;//Rinner,Router,alpha1->2 mirrored angle
  double Rradius;//mm receveir radius
  double Y0;//mm central point of one mirror composing the daisy
  double Beta=14.25/180.*pig;//tilt of one mirror composing the daisy
  double focal=6000.0; //mm
  double vi[3];//unit vector of incoming solar radiation
  double vr[3];//unit vector of reflected radiation
  double vn[3];//unit vector normal to mirror surface
  double P[3];//point of Dish's surface
  double P1[3];//impact point for -div on fluxmeter
  double P2[3];//impact point for +div on fluxmeter
  double PC[3];//impact point central on fluxmeter
  //FluxMeter
  int ng=201;//the FluxMeter is composed by ng x ng. pixel
  int nz=11;// for nz z levels
  double z1;//zmin  z-range analysed
  double z2;//zmax
  double zstep;
  double slg;//mm FluxMeter half-side
  double dg;//mm pixel side 
  double zg;//z plane where measuring flux
  double sd=0.00429;//<- Moon divergence instead of 0.00473; //rad, Sun divergence
  double g[201][201][3][11],gg[201][201];
  // g[ix][iy][0][iz]= pixel x coordinate
  // g[ix][iy][1][iz]= pixel y coordinate
  // g[ix][iy][2][iz]= concentration value
  double iad[360][90];//incidence angle distribution (phi 0-360,>>>> theta 0-90 <<<<)
  double iadmax;
  double SCR[91][100];
  Scalar scalar;
  int sumi;
  int njob=1;
  string pathbase="/home/marco/Workspace/Dish/";
  QString fileSCR="/home/marco/Workspace/RnearSpec/Almir/RvsTheta4.txt.SCR.txt";
  QString fileExpShape="/run/media/marco/d2quadra/VISdish/misura6nov2015/frames/results.txt";
  
  class ColorMap: public QwtLinearColorMap
{
public:
    ColorMap():
        QwtLinearColorMap(Qt::darkBlue, Qt::darkRed)
    {
        addColorStop(0.2, Qt::blue);
        addColorStop(0.4, Qt::cyan);
        addColorStop(0.6, Qt::yellow);
        addColorStop(0.8, Qt::red);
    }
};

  
SimulDish::SimulDish(QWidget *parent){
    setupUi(this); // this sets up GUI
 
    // signals/slots mechanism in action
    connect(pushButton_SCR,SIGNAL( clicked() ), this, SLOT( LoadSCR() ) );
    connect(pushButton_calculate, SIGNAL( clicked() ), this, SLOT( calculate() ) );
    connect(pushButton_reFlux, SIGNAL( clicked() ), this, SLOT( reFlux() ) );
    connect(pushButton_save, SIGNAL( clicked() ), this, SLOT( saveFlux() ) );
    connect(pushButton_load, SIGNAL( clicked() ), this, SLOT( loadFlux() ) );
    connect(pushButton_expFile, SIGNAL( clicked() ), this, SLOT( selectShape() ));
}

void SimulDish::closeEvent ( QCloseEvent *  ){
  qApp->quit();
}


void SimulDish::LoadSCR(){
 QString pathSCR="/home/marco/Workspace/RnearSpec/";
 fileSCR = QFileDialog::getOpenFileName(
        this,
        "Choose a Solar-radiation Conical Reflectance file",
        pathSCR,
        QString());
 lineEdit_SCR -> setText(fileSCR.section(pathSCR,1,1));
}


void SimulDish::selectShape(){
 QString pathExpShape="/run/media/marco/d2quadra/VISdish/";
 fileExpShape = QFileDialog::getOpenFileName(
        this,
        "Choose a Solar-radiation Conical Reflectance file",
        pathExpShape,
        QString());
 lineEdit_expFile -> setText(fileExpShape.section(pathExpShape,1,1));
}


int SimulDish::calculate(){
    double dr,da,ddr,rr,r,S,alpha,weight,phi,te,rh1,rh2,rmax,rmin,SW=0.,sai,sai2,cosinc,slopeDev;
    int Nrstep,Nastep,Nsstep,ixFMmin=ng-1,ixFMmax=0,iyFMmin=ng-1,iyFMmax=0,iphi,itheta,ithetainc;
    QString line;// = stream.readLine();

//parameter setting....
  focal = dSB_Rc -> value();
  focal=focal*1000;// focal in mm
  Qt::CheckState state1;
  nz = spinBox_Nz -> value();// z levels
  z1 = dSB_zmin -> value();
  z1=z1*1000;
  z2 = dSB_zmax -> value();
  z2=z2*1000;
  zstep=(z2-z1)/double(nz-1);
  slg = dSB_fluxmeter -> value();
  slg=slg/2.*1000.;
  dg=slg/double((ng-1)/2);
  Rradius = dSB_Dreceiver -> value();
  Rradius=Rradius/2.;
  S1P2 = comboBox_shape -> currentIndex();
  S1P2++;
  printf("S1P2= %d\n",S1P2);
  if(S1P2<3){
      state1 = checkBox_daisy -> checkState();
      if( state1 == Qt::Checked )
          daisy=1;
      else
          daisy=0;
      r1 = dSB_Finner -> value();
      r1=r1/2.*1000.;//inner radius in mm
      r2 = dSB_Fouter -> value();
      r2=r2/2.*1000.;//outer radius in mm
      if(r1>=r2){
          QMessageBox::critical(0, "Error", "Please set outer greater than inner!");
          return(1);
      }
      Y0=-1.425*r2;//a bit more than sqrt(2)=1.414 of r2
      double alphaDark = dSB_Dalpha -> value();
      alpha1=0.0;
      alpha2=(360.-alphaDark)*pig/180.;
      Beta = dSB_Beta -> value();
      Beta=Beta*pig/180;//Beta in rad
      ds = dSB_step -> value();// mm discretization step
      if(daisy!=1){
          vi[0]=0.0;
          vi[1]=0.0;
          vi[2]=-1.0;
      }
      else{
          vi[0]=0.0;
          vi[1]=sin(Beta);
          vi[2]=-cos(Beta);
      }
      Nrstep=int((r2-r1)/ds+0.5);
      dr=(r2-r1)/Nrstep;
  }
  else{//S1P2=3 = expShape
      daisy=0;
      vi[0]=0.;
      vi[1]=0.;
      vi[2]=-1.;
      Nrstep=1;
// diametroMaggiore = 11554 mm misurato il 16 Novembre 2015 da Pino
//      ds=11554./double(980);//reduction=2
      ds=11554./double(490);//reduction=4
//      ds=11554./double(196);//reduction=10
      dSB_step -> setValue(ds);
  }
//FluxMeter initialization
  for(int igz=0;igz<nz;igz++){ 
    for(int igx=0; igx<ng; igx++){
      for(int igy=0; igy<ng; igy++){
        g[igx][igy][0][igz]=-slg+dg*double(igx);//warning g has his own reference
        g[igx][igy][1][igz]=-slg+dg*double(igy);
        g[igx][igy][2][igz]=0.0;
      }
    }
  } 

// incidence angle distribution initialization
  for(int ii=0;ii<360;ii++){ 
    for(int jj=0; jj<90; jj++){
     iad[ii][jj]=0.0;
    }
  }
  
// load Solar-radiation Conical Reflectance
  QFile fSCR(fileSCR);
  int nphi;
  if (fSCR.open(QIODevice::ReadOnly | QIODevice::Text)){
    QTextStream stream ( &fSCR );
    line = fSCR.readLine();
    nphi=line.count("\t");
    printf("SCR=%s consists of %d columns\n",(fileSCR.toStdString()).c_str(),nphi);
    for(int i=0;i<=90;i++){
     for(int j=0;j<=nphi;j++){
      stream>>SCR[i][j];
      //cout<<SCR[i][j]<<"\t";
     }
     //cout<<"\n\n";
    }
    fSCR.close();
  }
  else
   nphi=1;//perfect solar-mirror
   
  state1 = checkBox_Rh -> checkState();
  if( state1 == Qt::Checked ) nphi=1;//hemispherical reflectance

  iadmax=0;
  sumi=0;
  sai=0.;
  sai2=0;
  printf("Computing with S1P2=%d....\n",S1P2);
  fflush(stdout);

  //open here the fileExpShape to avoid compiling errors
  QFile fileES(fileExpShape);
  if (!fileES.open(QIODevice::ReadOnly | QIODevice::Text)){
     printf("Attenzione: errore di apertura del file %s\n",(fileExpShape.toStdString()).c_str());
     return 17;
  }
  QTextStream stream2 ( &fileES );
  
  for(int i=0;i<Nrstep;i++){
    if(S1P2<3){
        printf("\rjob done(%%)= %d",int(100.*double(i*i)/double((Nrstep-1)*(Nrstep-1))+0.5));
        progressBar -> setValue(int(100.*double(i*i)/double((Nrstep-1)*(Nrstep-1))+0.5));
        fflush(stdout);
        r=r1+i*dr;
        Nastep=int(2.*(alpha2-alpha1)*r/ds+0.5);
        if(Nastep < 10) Nastep=10;
        da=(alpha2-alpha1)/Nastep;
        Nsstep=int(dr+0.5);//discretization for element-surface computing (about 1mm)
        ddr=dr/Nsstep;
        S=0.;
        for(int k=0;k<Nsstep;k++){//area of the 3D surface-element
            rr=r+k*ddr;
            S=S+da*rr*ddr*sqrt(1.+Dr(r)*Dr(r));
        }
        //    printf("r=%f Nastep=%d da=%f S=%f\n",r,Nastep,da,S);
    }
    else{
        Nastep=0;
        while ( !stream2.atEnd() ) {
            Nastep++;
            line = stream2.readLine( );
        }
        printf("The selected experimentalShape file consists of Nastep = %d data\n",Nastep);
        fileES.seek(0);//rewind stream2
    }
    for(int j=0;j<Nastep;j++){
        cout << j << " of "<<Nastep<<"\n";
        fflush(stdout);
      if(S1P2<3){
          //      printf("r=%f done(%%)=%f\n",r,100.*float(j)/float(Nastep));
          alpha=alpha1+j*da;
          cosinc=reflection(r,alpha);
          ithetainc=int(acos(cosinc)*180./pig+0.5);
          weight=S*cosinc;//each element is weighted for S*cos(theta_inc)
          SW=SW+weight;
          P[0]=r*cos(alpha);
          P[1]=r*sin(alpha);
          if(S1P2==1)
              P[2]=2.*focal-sqrt(4.*focal*focal-r*r);
          else
              P[2]=r*r/4.0/focal;
          if(daisy==1) TrasformXYZ(Y0,Beta);
      }
      else{
         progressBar -> setValue(int(100.*double(j)/double(Nastep)+0.5));
         int jskip=-1;//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 3395 mm è il raggio del bordo pannelli interni
         do{
            stream2 >> P[0] >> P[1] >> P[2] >> vn[0] >> vn[1] >> vn[2] >> slopeDev >>slopeDev >> slopeDev;
            jskip++;
         } while(sqrt(P[0]*P[0]+P[1]*P[1])>3395. && j+jskip<Nastep);//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< per limitare ai soli pannelli interni
         j=j+jskip;
         cosinc=unitVectorReflectance();
         ithetainc=int(acos(cosinc)*180./pig+0.5);
         weight=ds*ds;//the pixel surface is parallel to XY plane so cosheta_inc=1
         SW=SW+weight;
//         printf("P: %f %f %f\n",P[0],P[1],P[2]);
//         printf("vi: %f %f %f\n",vi[0],vi[1],vi[2]);
//         printf("vn: %f %f %f\n",vn[0],vn[1],vn[2]);
//         printf("vr: %f %f %f\n",vr[0],vr[1],vr[2]);
//         sleep(10);
      }

      phi=atan2(vr[1],vr[0]);// polar coordinate angles
      te=acos(vr[2]);
      if(phi > 0)
       iphi=int(phi*180./pig);
      else
       iphi=int((phi+2.*pig)*180./pig);
      itheta=int(te*180./pig);
      iad[iphi][itheta]++;
      if(iad[iphi][itheta] > iadmax) iadmax=iad[iphi][itheta];
      sumi++;
      sai=sai+te;
      sai2=sai2+te*te;
      
//      printf("vr:%f, %f, %f phi=%f te=%f\n",vr[0],vr[1],vr[2],phi*180./pig,te*180./pig);
      for(int igz=0;igz<nz;igz++){
        zg=z1+igz*zstep;
        double ssd,dSCR,ee0,A0=0.,A,rmin0=0.,rmax0=0.;
        for(int iPHI=1;iPHI<=nphi;iPHI++){
          if(nphi==1){
           ssd=sd;
           dSCR=SCR[ithetainc][0];//=Hemispherical reflectance
          }
          else{
           ssd=double(iPHI)/1000.;
           if(iPHI==1)
            dSCR=SCR[ithetainc][iPHI];
           else
            dSCR=SCR[ithetainc][iPHI]-SCR[ithetainc][iPHI-1];
          }
          rh1=(zg-P[2])/cos(te+ssd);// path length to FluxMeter for +div
          P1[0]=rh1*sin(te+ssd)*cos(phi)+P[0];
          P1[1]=rh1*sin(te+ssd)*sin(phi)+P[1];
          P1[2]=zg;
          rh2=(zg-P[2])/cos(te-ssd);// path length to FluxMeter for -div
          P2[0]=rh2*sin(te-ssd)*cos(phi)+P[0];
          P2[1]=rh2*sin(te-ssd)*sin(phi)+P[1];
          P2[2]=zg;
          PC[0]=(P1[0]+P2[0])/2.;
          PC[1]=(P1[1]+P2[1])/2.;
          PC[2]=zg;
          rmax=0.5*sqrt((P1[0]-P2[0])*(P1[0]-P2[0])+(P1[1]-P2[1])*(P1[1]-P2[1]));
          rmin=ssd*sqrt((P[0]-PC[0])*(P[0]-PC[0])+
                      (P[1]-PC[1])*(P[1]-PC[1])+
                      (P[2]-PC[2])*(P[2]-PC[2]));
          A=pig*rmin*rmin;          
    //       printf(" P: %f, %f %f\n P1: %f, %f %f\n P2: %f, %f %f\n",P[0],P[1],P[2],
    //              P1[0],P1[1],P1[2],P2[0],P2[1],P2[2]);
    //      printf("rmin=%f rmax=%f\n",rmin,rmax);
          if(rmax<rmin) {
            printf("rmax<rmin!\n");
    //        return(0);
          }
          int igxMIN=int((PC[0]-rmax*1.05+slg)/dg);
          int igxMAX=int((PC[0]+rmax*1.05+slg)/dg+0.5);
          int igyMIN=int((PC[1]-rmax*1.05+slg)/dg);
          int igyMAX=int((PC[1]+rmax*1.05+slg)/dg+0.5);
          ixFMmin=min(ixFMmin,igxMIN);
          ixFMmax=max(ixFMmax,igxMAX);
          iyFMmin=min(iyFMmin,igyMIN);
          iyFMmax=max(iyFMmax,igyMAX);
          if(igxMIN < 0) igxMIN=0;
          if(igxMIN >= ng) igxMIN=ng-1;
          if(igxMAX < 0)  igxMAX=0;
          if(igxMAX >= ng) igxMAX=ng-1;
          if(igyMIN < 0) igyMIN=0;
          if(igyMIN >= ng) igyMIN=ng-1;
          if(igyMAX < 0)  igyMAX=0;
          if(igyMAX >= ng) igyMAX=ng-1;
          double xc,yc,mv,xr,yr,ee,phic;
          for(int igx=igxMIN; igx<=igxMAX; igx++){
            for(int igy=igyMIN; igy<=igyMAX; igy++){
              xc=g[igx][igy][0][igz]-PC[0];
              yc=g[igx][igy][1][igz]-PC[1];
              mv=sqrt(xc*xc+yc*yc);
              phic=atan2(yc,xc);
              xr=mv*cos(phic-phi);
              yr=mv*sin(phic-phi);
              ee=(xr/rmax)*(xr/rmax)+(yr/rmin)*(yr/rmin);
              if(iPHI>1)
               ee0=(xr/rmax0)*(xr/rmax0)+(yr/rmin0)*(yr/rmin0);
              else
               ee0=2.;
              if(ee0 >=1.0 && ee <= 1.0){
                g[igx][igy][2][igz]=g[igx][igy][2][igz]+weight/(A-A0)*cos(te)*dSCR;
    //             vp(1)=g(igx,igy,1)-pe(1)
    //             vp(2)=g(igx,igy,2)-pe(2)
    //             vp(3)=focale-pe(3)
    //             m2=vp(1)*vp(1)+vp(2)*vp(2)+vp(3)*vp(3)
    //             azi=atan2(vp(2),vp(1))
    //             zen=acos(vp(3)/sqrt(m2))
    //             i=int((azi+pig)/bai)+1
    //             j=int(zen/bai)+1
    //             dai(i,j)=dai(i,j)+pai*cos(te)*peso
              }
            }
          }
          rmin0=rmin;
          rmax0=rmax;
          A0=A;
        }
      }
    }
  }

//daisy flux overlapping 
  if(daisy==1){
    SW=SW*4.;
    int ngm=(ng-1)/2;
    for(int k=0;k<nz;k++){
      for(int i=0; i<=ngm; i++){
        for(int j=0; j<=ngm; j++){
          gg[i][j]=g[i][j][2][k]+g[ng-1-j][i][2][k]+
                        g[j][ng-1-i][2][k]+g[ng-1-i][ng-1-j][2][k];             
        }
      }
      for(int i=0; i<=ngm; i++){
        for(int j=0; j<=ngm; j++){
          g[i][j][2][k]=gg[i][j];
          g[ng-1-i][j][2][k]=gg[i][j];
          g[i][ng-1-j][2][k]=gg[i][j];
          g[ng-1-j][ng-1-i][2][k]=gg[i][j];
        }
      }
    }
  }

// Flux analysis
  int NpxReceiver,izo=0;
  double Cmax,CmeanCaptured,SumTot,SumIntercepted,dist,FMZ[100][5],RatioCaptured,RCmax=0.;
  for(int k=0;k<nz;k++){
    NpxReceiver=0;
    Cmax=0.0;
    CmeanCaptured=0.0;
    SumTot=0.0;
    SumIntercepted=0.0;
    for(int i=0; i<ng; i++){
      for(int j=0; j<ng; j++){
        SumTot=SumTot+g[i][j][2][k];
        dist=sqrt(g[i][j][0][k]*g[i][j][0][k]+g[i][j][1][k]*g[i][j][1][k]);
        if(dist<=Rradius){
          SumIntercepted=SumIntercepted+g[i][j][2][k];
          NpxReceiver++;
          Cmax=max(Cmax,g[i][j][2][k]);
        }
      }
    }
    CmeanCaptured=SumIntercepted/double(NpxReceiver);
    RatioCaptured=SumIntercepted/SumTot*100.;
    FMZ[k][0]=z1+k*zstep;
    FMZ[k][1]=RatioCaptured;
    FMZ[k][2]=Cmax;
    FMZ[k][3]=CmeanCaptured;
    FMZ[k][4]=CmeanCaptured*pig*Rradius*Rradius/1.0e+06;//(suns*m^2)
    if(RatioCaptured > RCmax){
     RCmax=RatioCaptured;
     izo=k;
     zg=FMZ[k][0];
    }
  }
  
// normalization of incidence angle distribution and statistics
  for(int i=0;i<360;i++){
   for(int j=0;j<90;j++){
    iad[i][j]=iad[i][j]/iadmax;
   }
  }
 
// check njob to be used for file saving 
  char fileFluxMap[200],fileiad[200];
  int control=0,ierr;
  while(control==0){
   ierr = sprintf(fileFluxMap,"%sfluxMap_job_%d.jpeg",pathbase.c_str(),njob);
   if(ierr<0) printf("An error occurred at kine 419\n");
   QFile file(fileFluxMap);
   if (!file.exists()){
    control=1;
   }
   else{
    njob++;
    lineEdit_Njob-> setText(QString::number(njob));
   }
  }
  ierr = sprintf(fileiad,"%sIAD_job_%d.jpeg",pathbase.c_str(),njob);//file incidence angle distribution
  
//message
  QString msg;
  if(daisy==1 && S1P2==1)
    printf("\nDaisy dish composed by 4 spherical mirrors\n");
  else if(daisy==1 && S1P2==2)
    printf("\nDaisy dish composed by 4 paraboloidal mirrors\n");
  else if(daisy!=1 && S1P2==1)
    printf("\nSpherical dish\n");
  else if(daisy!=1 && S1P2==2)
    printf("\nParaboloidal dish\n");
  if(daisy==1)
    printf("Y0= %f Beta(deg)= %f\n",Y0,Beta/pig*180.);
  printf("FluxMeter: side=%f mm  pixel=%f mm in z= %f mm\n",2.*slg,dg,zg);
  printf("Nrstep=%d dr=%f\n",Nrstep,dr);
  printf("Effective surface (m^2)= %f\n",SW/1.E6);
  if(ixFMmin < 0) {
    printf("Please enlarge FluxMeter to xMIN=%f\n",-slg+dg*(ixFMmin));
    msg="Please enlarge FluxMeter to xMIN= "+QString::number((-slg+dg*(ixFMmin))/1000.);
    QMessageBox::critical(0, "Warning", msg);
  }
  if(ixFMmax >= ng){
    printf("Please enlarge FluxMeter to xMAX=%f\n",slg+dg*(ixFMmax-ng));
    msg="Please enlarge FluxMeter to xMAX= "+QString::number((slg+dg*(ixFMmax-ng))/1000.);
    QMessageBox::critical(0, "Warning", msg);
  }
  if(iyFMmin < 0) {
    printf("Please enlarge FluxMeter to yMIN=%f\n",-slg+dg*(iyFMmin));
    msg="Please enlarge FluxMeter to yMIN= "+QString::number((-slg+dg*(iyFMmin))/1000.);
    QMessageBox::critical(0, "Warning", msg);
  }
  if(iyFMmax >= ng){
    printf("Please enlarge FluxMeter to yMAX=%f\n",slg+dg*(iyFMmax-ng));
    msg="Please enlarge FluxMeter to yMAX= "+QString::number((slg+dg*(iyFMmax-ng))/1000.);
    QMessageBox::critical(0, "Warning", msg);
  }
  printf("ixFMmax=%d ixFMmin=%d iyFMmax=%d iyFMmin=%d\n",ixFMmax,ixFMmin,iyFMmax,iyFMmin);
  int isopti=max(ixFMmax-ixFMmin,iyFMmax-iyFMmin);
  if(double(isopti)/double(ng) < 0.9){
    printf("Please reduce FluxMeter side to %f\n",dg*double(isopti)/1000.);
    msg="Please reduce FluxMeter side to "+QString::number(dg*double(isopti)/1000.)+" m";
    QMessageBox::critical(0, "Warning", msg);
  }
  printf("Optimal fluxmeter-side= %f m\n",dg*double(isopti)/1000.);
  printf("Receiver: radius= %f z_optimal= %f\n",Rradius,zg);
  printf("RatioCaptured(%%)= %f\nCmax= %f\nCmean_captured= %f\n",
         FMZ[izo][1],FMZ[izo][2],FMZ[izo][3]);
  lineEdit_Njob-> setText(QString::number(njob));
  lineEdit_A -> setText(QString::number(SW/1.E6));
  lineEdit_zOpt -> setText(QString::number(zg/1000.));
  lineEdit_RC -> setText(QString::number(FMZ[izo][1]));
  lineEdit_Cmax -> setText(QString::number(FMZ[izo][2]));
  lineEdit_Cmean -> setText(QString::number(FMZ[izo][3]));
  lineEdit_yield -> setText(QString::number(FMZ[izo][4]));
  lineEdit_Ceffi -> setText(QString::number(FMZ[izo][4]/(SW/1.E6)));
  lineEdit_thetaInc -> setText(QString::number(sai/double(sumi)*180./pig)+QChar(0xb1)+
           QString::number(sqrt(sai2/sumi-(sai/double(sumi))*(sai/double(sumi)))*180./pig));
  spinBox_NzF -> setValue(izo);
  
  //PLOTs
  double xPlot[nz],yPlot[nz],yCmean[nz],yCmax[nz];
  for(int k=0;k < nz;k++){
    xPlot[k]=z1+k*zstep;
    yPlot[k]=FMZ[k][1];
    yCmax[k]=FMZ[k][2];
    yCmean[k]=FMZ[k][3];
  }
  QwtPlot *RCPlot=new QwtPlot();
  QwtPlotCurve *curve1=new QwtPlotCurve("Curve 1");
  RCPlot -> setTitle("RC - #job= "+QString::number(njob));
  RCPlot -> setAxisTitle(0,"RC (%)");
  RCPlot -> setAxisTitle(2,"z (mm)");
  curve1->setSamples(xPlot, yPlot, nz);
  curve1->attach(RCPlot);
  // Make the grid on
  QwtPlotGrid *grid = new QwtPlotGrid();
  grid->setPen(QPen(Qt::gray, 0.0, Qt::DotLine));
  grid->enableX(true);
  grid->enableXMin(true);
  grid->enableY(true);
  grid->enableYMin(true);
  grid->attach(RCPlot);
  RCPlot->replot();
  RCPlot->show();
  
  QwtPlot *CPlot=new QwtPlot();
  QwtPlotCurve *curve2=new QwtPlotCurve("Curve 2");
  CPlot -> setTitle("Concentration - #job= "+QString::number(njob));
  CPlot -> setAxisTitle(0,"C (Suns)");
  CPlot -> setAxisTitle(2,"z (mm)");
  curve2->setSamples(xPlot, yCmax, nz);
  curve2->setPen(QPen(Qt::red,2,Qt::SolidLine));
  curve2->attach(CPlot);
  QwtPlotCurve *curve3=new QwtPlotCurve("Curve 3");
  curve3->setSamples(xPlot, yCmean, nz);
  curve3->setPen(QPen(Qt::blue,2,Qt::DotLine));
  curve3->attach(CPlot);
   
  // Make the grid on
  QwtPlotGrid *grid2 = new QwtPlotGrid();
  grid2->setPen(QPen(Qt::gray, 0.0, Qt::DotLine));
  grid2->enableX(true);
  grid2->enableXMin(true);
  grid2->enableY(true);
  grid2->enableYMin(true);
  grid2->attach(CPlot);
  CPlot->replot();
  CPlot->show();
  
  // contour map of flux
  Mat smap(ng,ng,CV_8UC3);// = cvCreateImage( cvSize( ng, ng ), IPL_DEPTH_8U, 3 );
  double flux;
  for(int irs=0;irs<ng;irs++){
    for(int jrs=0;jrs<ng;jrs++){
      flux=g[jrs][ng-1-irs][2][izo]/FMZ[izo][2];
      setBGR(flux);
      smap.at<Vec3s>(irs,jrs)[0]=scalar[0];
      smap.at<Vec3s>(irs,jrs)[1]=scalar[1];
      smap.at<Vec3s>(irs,jrs)[2]=scalar[2];
    }
  }
  Point2i Pt1;
  Pt1.x=int(double(ng-1)/2.);
  Pt1.y=int(double(ng-1)/2.);
  circle(smap,Pt1,int(Rradius/dg+0.5),Scalar(255,255,255),1,8,0);
  imwrite(fileFluxMap, smap );
  namedWindow("flux map");
  moveWindow("flux map", 510, 650);
  imshow("flux map", smap );
  smap.release();
  
//   QwtPlot *fluxPlot=new QwtPlot();
//   QwtPlotCurve *curve4=new QwtPlotCurve("Curve 4");
//   fluxPlot -> setTitle("Flux - #job= "+QString::number(njob));
//   fluxPlot -> setAxisTitle(0,"x (mm)");
//   fluxPlot -> setAxisTitle(2,"y (mm)");
//   curve4->setSamples(xPlot, yCmax, nz);
//   curve4->setSymbol(new QwtSymbol(QwtSymbol::Ellipse, Qt::yellow,
//         QPen(Qt::yellow), QSize(5, 5) ) );
//   curve4->setStyle(QwtPlotCurve::NoCurve);
//   curve4->attach(fluxPlot);
//   fluxPlot->replot();
//   fluxPlot->show();
 
//   QVector<double> gv;
//   for(int i=0; i<ng; i++){
//     for(int j=0; j<ng; j++){
//        gv << g[j][i][2][izo];
//     }
//   }
//   QwtPlot *CfluxPlot=new QwtPlot();
//   setCentralWidget(CfluxPlot);
//   d_spectrogram = new QwtPlotSpectrogram();
//   d_spectrogram->setRenderThreadCount(0); // use system specific thread count
//   d_spectrogram->setColorMap( new ColorMap() );
//   QwtMatrixRasterData *RasterData = new QwtMatrixRasterData();
//   RasterData->setValueMatrix(gv,ng);
//   d_spectrogram->setData(RasterData);
//   d_spectrogram->attach(CfluxPlot);
// 
//   const QwtInterval zInterval = d_spectrogram->data()->interval( Qt::ZAxis );
//   // A color bar on the right axis
//   QwtScaleWidget *rightAxis = axisWidget(QwtPlot::yRight);
//   rightAxis->setColorBarEnabled(true);
//   rightAxis->setColorBarWidth(40);
//   rightAxis->setColorMap(zInterval, new ColorMap() );
// 
//   setAxisScale(QwtPlot::yRight, zInterval.minValue(), zInterval.maxValue() );
//   enableAxis(QwtPlot::yRight);
// 
//   plotLayout()->setAlignCanvasToScales(true);
// 
//   setAxisScale(QwtPlot::xBottom, -slg, slg);
//   setAxisMaxMinor(QwtPlot::xBottom, 0);
//   setAxisScale(QwtPlot::yLeft, -slg, slg);
//   setAxisMaxMinor(QwtPlot::yLeft, 0);
// 
//   QwtPlotMagnifier *magnifier = new QwtPlotMagnifier( canvas() );
//   magnifier->setAxisEnabled( QwtPlot::yRight, false);
// 
//   QwtPlotPanner *panner = new QwtPlotPanner( canvas() );
//   panner->setAxisEnabled( QwtPlot::yRight, false);
// 
//   canvas()->setBorderRadius( 10 );
//   

// contour map of incidence angle distribution
  //****************  creation of the incidence angle map
  Mat smap2( 360, 90 , CV_8UC3 );
  //step      = smap2->widthStep;
  //chan  = smap2->nChannels;
  //uchar *datasmap2    = (uchar *)smap2->imageData;
  double iav;
  for(int irs=0;irs<360;irs++){
    for(int jrs=0;jrs<90;jrs++){
      iav=iad[irs][jrs];
      setBGR(iav);
      smap2.at<Vec3s>(irs,jrs)[0]=scalar[0];
      smap2.at<Vec3s>(irs,jrs)[1]=scalar[1];
      smap2.at<Vec3s>(irs,jrs)[2]=scalar[2];
    }
  }
  imwrite(fileiad, smap2 );
  namedWindow("flux map2");
  moveWindow("flux map2", 510, 650);
  imshow("flux map2", smap2 );
  smap2.release();
  
  // save the flux map at the optimal z
  ierr = sprintf(fileFluxMap,"%sfluxMap_job_%d.txt",pathbase.c_str(),njob);
  ofstream ffm;
  ffm.open(fileFluxMap);
  if(daisy==1 && S1P2==1){
    ffm << "Daisy dish composed by 4 spherical mirrors" <<"\n";
    ffm << "Focal length (mm) " << focal <<"\n";
  }
  else if(daisy==1 && S1P2==2){
    ffm << "Daisy dish composed by 4 paraboloidal mirrors" <<"\n";
    ffm << "Focal length (mm) " << focal <<"\n";
  }
  else if(daisy!=1 && S1P2==1){
    ffm << "Spherical dish" <<"\n";
    ffm << "Focal length (mm) " << focal <<"\n";
  }
  else if(daisy!=1 && S1P2==2){
   ffm << "Paraboloidal dish" <<"\n";
   ffm << "Focal length (mm) " << focal <<"\n";
  }
  if(daisy==1) {
   ffm << "Y0(mm) " << Y0 <<"\n";
   ffm << "Beta(deg) " << Beta/pig*180. <<"\n";
  }
  ffm << "radius_inner(mm) " << r1  << "\n";
  ffm << "radius_outer(mm) " << r2  << "\n";
  ffm << "alpha1(rad) " << alpha1  << "\n";
  ffm << "alpha2(rad) " << alpha2  << "\n";
  ffm << "Solar-radiation Conical Reflectance " << (fileSCR.toStdString()).c_str() << "\n";
  if(nphi == 1)
   ffm << "SCR = Rhemispherical in "<< sd <<" rad\n";
  else
   ffm << "SCR in (0,"<<nphi<<") mrad\n";
  ffm << "Receiver_radius(mm) " << Rradius << "\n";
  ffm << "Effective_surface(m^2) " << SW/1.E6 << "\n";
  ffm << "RatioCaptured(%) " << FMZ[izo][1] << "\n";
  ffm << "Cmax(suns) " << FMZ[izo][2] << "\n";
  ffm << "Cmean_captured(suns) " << FMZ[izo][3] << "\n";
  ffm << "Yield(suns*m^2)) " << FMZ[izo][4] << "\n";
  ffm << "Concentration_Efficiency " << FMZ[izo][4]/(SW/1.E6) << "\n";
  ffm << "FluxMeter_side(mm) " << 2.*slg << "\n";
  ffm << "FluxMeter_pixel_side(mm) " << dg << "\n";
  ffm << "FluxMeter_optimal_z(mm) " << zg << "\n";
  ffm << "Bin max (normalized to N. rays) of Incidence Angle Distribution " << iadmax/double(sumi) << "\n";
  ffm << "Incidence Angle (deg): mean = " << sai/double(sumi)*180./pig << "   sigma = "<< sqrt(sai2/sumi-(sai/double(sumi))*(sai/double(sumi)))*180./pig << "\n";
  ffm << "**** z, RatioCaptured, Cmax, CmeanCaptured *****\n";
  for(int k=0;k < nz;k++){
    ffm << z1+k*zstep <<"\t"<< FMZ[k][1]<<"\t"<< FMZ[k][2]<<"\t"<< FMZ[k][3] <<"\n";
  }
  ffm << "**** Concentration data at optimal z for 3D map drawing *****\n";
  for(int i=0; i<ng; i++){
    for(int j=0; j<ng; j++){
      if(j != ng-1)
        ffm << g[j][i][2][izo] << "\t";
      else
        ffm << g[j][i][2][izo] << "\n";
    }
  }
  ffm.close();
  
  // save the incidence angle distribution data
  ierr = sprintf(fileiad,"%sIAD_job_%d.txt",pathbase.c_str(),njob);
  ofstream fiad;
  fiad.open(fileiad);
  fiad << "**** Incidence Angle Distribution (theta x phi) *****\n";
  for(int i=0; i<360; i++){
    for(int j=0; j<90; j++){
      if(j != 89)
        fiad << iad[i][j] << "\t";
      else
        fiad << iad[i][j] << "\n";
    }
  }
  fiad.close();
  

  return(0);
}


void SimulDish::reFlux(){
 // Flux re-analysis
  int NpxReceiver;
  int k=spinBox_NzF -> value();
  double z=z1+k*zstep;
  lineEdit_zCurrent -> setText(QString::number(z));
  double Dx=dSB_Dx -> value();
  double Dy=dSB_Dy -> value();
  Rradius=dSB_Dreceiver -> value();
  Rradius=Rradius/2.;
  double A=lineEdit_A -> text().toDouble();
  double Cmax,CmeanCaptured,SumTot,SumIntercepted,dist,RatioCaptured;
  NpxReceiver=0;
  Cmax=0.0;
  CmeanCaptured=0.0;
  SumTot=0.0;
  SumIntercepted=0.0;
  for(int i=0; i<ng; i++){
    for(int j=0; j<ng; j++){
      SumTot=SumTot+g[i][j][2][k];
      dist=sqrt((g[i][j][0][k]-Dx)*(g[i][j][0][k]-Dx)+(g[i][j][1][k]-Dy)*(g[i][j][1][k]-Dy));
      if(dist<=Rradius){
        SumIntercepted=SumIntercepted+g[i][j][2][k];
        NpxReceiver++;
        Cmax=max(Cmax,g[i][j][2][k]);
      }
    }
  }
  CmeanCaptured=SumIntercepted/double(NpxReceiver);
  RatioCaptured=SumIntercepted/SumTot*100.;
  lineEdit_RC -> setText(QString::number(RatioCaptured));
  lineEdit_Cmax -> setText(QString::number(Cmax));
  lineEdit_Cmean -> setText(QString::number(CmeanCaptured));
  lineEdit_yield -> setText(QString::number(CmeanCaptured*pig*Rradius*Rradius/1.0e+06));
  lineEdit_Ceffi -> setText(QString::number(CmeanCaptured*pig*Rradius*Rradius/1.0e+06/A));
}


void SimulDish::saveFlux(){
 int ivalue;
 double value;
 Qt::CheckState state;
 QString filename,string;
 QString path="/home/marco/Workspace/Dish";
 filename=QFileDialog::getSaveFileName(
          this,
          "Filename to save",
          path,
          "simulDish Flux (*.sdf)");
 printf("filename = %s\n",(filename.toStdString()).c_str());
 QFile fileS(filename);
 if(!fileS.open(QIODevice::WriteOnly | QIODevice::Text))
  return;
 QTextStream out(&fileS);
 value=dSB_Rc -> value();
 out << "focalLength(m)\t" << value << "\n";
 value=dSB_Fouter -> value();
 out << "Douter(m)\t"<< value << "\n";
 value=dSB_Finner -> value();
 out << "Dinner(m)\t"<< value << "\n";
 value=dSB_Dalpha -> value();
 out << "darkSlice(deg)\t"<< value << "\n";
 value=dSB_Beta -> value();
 out << "Beta(deg)\t"<< value << "\n";
 value=dSB_Dreceiver -> value();
 out << "Dreceiver(mm)\t"<< value << "\n";
 value=dSB_zmin -> value();
 out << "zmin(m)\t"<< value << "\n";
 value=dSB_zmax -> value();
 out << "zmax(m)\t"<< value << "\n";
 value=dSB_fluxmeter -> value();
 out << "Side(m)\t"<< value << "\n";
 ivalue=spinBox_Nz -> value();
 out << "NzLevels\t" << ivalue << "\n";
 ivalue=comboBox_shape -> currentIndex();
 ivalue++;
 out << "S1P2\t" << ivalue << "\n";
 state=checkBox_daisy -> checkState();
 if(state==Qt::Unchecked)
   ivalue=0;
 else
   ivalue=1;
 out << "Daisy\t" << ivalue << "\n";
 state=checkBox_Rh -> checkState();
 if(state==Qt::Unchecked)
   ivalue=0;
 else
   ivalue=1;
 out << "SCR=Rh\t" << ivalue << "\n";
 out << "fileSCR\t" << fileSCR << "\n";
 value=dSB_step -> value();
 out << "DisStep(mm)\t"<< value << "\n";
 string=lineEdit_A -> text();
 value=string.toDouble();
 out << "Aeffective(m^2)\t" << value << "\n";
 string=lineEdit_zOpt -> text();
 value=string.toDouble();
 out << "zOptimal(m)\t" << value << "\n";
 ivalue=spinBox_NzF -> value();
 out << "zLevelDisplayed\t" << ivalue << "\n";
 out << "zLevelTot\t" << nz << "\n";
 out << "ng\t" << ng << "\n";
 for(int igz=0;igz<nz;igz++){ 
    for(int igx=0; igx<ng; igx++){
      for(int igy=0; igy<ng; igy++){
       out << g[igx][igy][0][igz] <<"\t"<<g[igx][igy][1][igz]<<"\t"<<g[igx][igy][2][igz];
       if(igy<ng-1)
        out << "\t";
       else
        out << "\n";
      }
    }
 }
 fileS.close();
}


void SimulDish::loadFlux(){
 int ivalue;
 double value;
 QString filename,string;
 QString path="/home/marco/Workspace/Dish";
 filename = QFileDialog::getOpenFileName(
        this,
        "Choose a SimulDish Flux file",
        path,
        "simulDish Flux (*.sdf)");
 QFile fileL(filename);
 if (!fileL.open(QIODevice::ReadOnly | QIODevice::Text)) 
  return;
 QTextStream stream (&fileL);
 stream >> string >> value;
 dSB_Rc -> setValue(value);
 stream >> string >> value;
 dSB_Fouter -> setValue(value);
 stream >> string >> value;
 dSB_Finner -> setValue(value);
 stream >> string >> value;
 dSB_Dalpha -> setValue(value);
 stream >> string >> value;
 dSB_Beta -> setValue(value);
 stream >> string >> value;
 dSB_Dreceiver -> setValue(value);
 stream >> string >> value;
 dSB_zmin -> setValue(value);
 z1=value;
 stream >> string >> value;
 dSB_zmax -> setValue(value);
 z2=value;
 stream >> string >> value;
 dSB_fluxmeter -> setValue(value);
 stream >> string >> ivalue;
 spinBox_Nz -> setValue(value);
 zstep=(z2-z1)/(ivalue-1);
 stream >> string >> ivalue;
 comboBox_shape -> setCurrentIndex(ivalue-1);
 stream >> string >> ivalue;
 if(ivalue==0)
  checkBox_daisy -> setCheckState(Qt::Unchecked);
 else
  checkBox_daisy -> setCheckState(Qt::Checked);
 stream >> string >> ivalue;
 if(ivalue==0)
  checkBox_Rh -> setCheckState(Qt::Unchecked);
 else
  checkBox_Rh -> setCheckState(Qt::Checked);
 stream >> string >> fileSCR;
 stream >> string >> value;
 dSB_step -> setValue(value);
 stream >> string >> value;
 lineEdit_A ->setText(QString::number(value));
 stream >> string >> value;
 lineEdit_zOpt ->setText(QString::number(value));
 stream >> string >> ivalue;
 spinBox_NzF -> setValue(value);
 stream >> string >> nz;
 stream >> string >> ng;
 for(int igz=0;igz<nz;igz++){ 
    for(int igx=0; igx<ng; igx++){
      for(int igy=0; igy<ng; igy++){
       stream >> g[igx][igy][0][igz] >> g[igx][igy][1][igz] >> g[igx][igy][2][igz];
      }
    }
 }
 fileL.close();
}

double Dr(double r){
 double ratio=0.5*r/focal;
 if(S1P2==1)
   ratio=ratio/sqrt(1.-ratio*ratio);//Error here: was * insteat of / !!!
 return(ratio);
}


double reflection(double r, double alpha){
  double vn[3],mvn=0.,cosinc=0.;
  vn[0]=-Dr(r)*cos(alpha);
  vn[1]=-Dr(r)*sin(alpha);
  vn[2]=1.0;
  for(int i=0;i<3;i++){
    mvn=mvn+vn[i]*vn[i];
    cosinc=cosinc+vi[i]*vn[i];
  }
  mvn=sqrt(mvn);
  cosinc=-cosinc/mvn;
  for(int i=0;i<3;i++){
   vr[i]=vi[i]+2.0*cosinc*vn[i]/mvn;
  }
//   printf("reflection:\n vi:%f, %f, %f\n vr:%f, %f, %f\n",
//          vi[0],vi[1],vi[2],vr[0],vr[1],vr[2]);
  return(cosinc);
}


double unitVectorReflectance(){
    double cosinc=0.,mvn=0.;
    for(int i=0;i<3;i++)
      cosinc=cosinc+vi[i]*vn[i];
    cosinc=-cosinc;
    for(int i=0;i<3;i++){
     vr[i]=vi[i]+2.0*cosinc*vn[i];
     mvn=mvn+vr[i]*vr[i];
    }
    for(int i=0;i<3;i++)
        vr[i]=vr[i]/sqrt(mvn);
  //   printf("reflection:\n vr:%f, %f, %f\n",
  //          vr[0],vr[1],vr[2]);
    return(cosinc);
}


void TrasformXYZ(double yt, double bet){
  double Pt[3],vrt[3];
  Pt[0]=P[0];
  Pt[1]=P[1]*cos(bet)+P[2]*sin(bet)+yt;
  Pt[2]=-P[1]*sin(bet)+P[2]*cos(bet);
  vrt[0]=vr[0];
  vrt[1]=vr[1]*cos(bet)+vr[2]*sin(bet);
  vrt[2]=-vr[1]*sin(bet)+vr[2]*cos(bet);
  for(int i=0;i<3;i++){
   P[i]=Pt[i];
   vr[i]=vrt[i];
  }
}


void setBGR(double flux){
    int Blue=0;
    int Green=0;
    int Red=0;
    if(flux>0.0 && flux< 0.25){
     Blue=255;
     Green=int(255.*4.*flux+0.5);
     Red=0;
    }
    else if(flux>= 0.25 && flux< 0.50){
     Blue=int(255.*(1.-4.*(flux-0.25)));
     Green=255;
     Red=0;
    }
    else if(flux>= 0.50 && flux< 0.75){
     Blue=0;
     Green=255;
     Red=int(255.*4.*(flux-0.50));
    }
    else if(flux>= 0.75){
     Blue=0;
     Green=int(255.*(1.-4.*(flux-0.75)));
     Red=255;
    }
    if(Green<0) Green=0;
    scalar.val[0]=Blue;
    scalar.val[1]=Green;
    scalar.val[2]=Red;
}
