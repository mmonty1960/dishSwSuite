/*Author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it

ImageWorkBench processes the images acquired by kFluxMapper
returning the meta-image (double float) mapping the experimental flux


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
#include "iwb.h"
#include <QtGui>
#include <QFile>
#include <QFileDialog>
#include <fstream>
#include <iostream>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/calib3d/calib3d.hpp"

using namespace std;
using namespace cv;

/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia
Roma - Italy*/

//global variables *********************************
//QString dir="/run/media/marco/d2quadra/FluxMapper";//"/run/media/marco/d2quadra/FluxMapper/moonTracking20mag2016/z41cm"
QString dir="/home/marco/Workspace/FluxMapper";
QString filename;
int mpxX1,mpxY1,mpxX2,mpxY2,drag;
int iMIN,iMAX,jMIN,jMAX;
double SumFlux[100][4];
double pxImDim=0.4;//mm
Mat img;
Mat Cimg(2050,2448,CV_64FC1);//composite image
Mat imgSt;

//invoked functions ********************************
void mouseHandler(int event, int x, int y, int flags, void* param);


ImageWorkBench::ImageWorkBench(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::ImageWorkBench)
{
    ui->setupUi(this);

    ui->lineEdit_dir -> setText(dir);

    connect( ui->pB_setDir,    SIGNAL( clicked() ),               this, SLOT( selectDir() ) );
    connect( ui->pB_compute,   SIGNAL( clicked() ),               this, SLOT( compute() ) );
    connect( ui->pB_Load,      SIGNAL( clicked() ),               this, SLOT( selectImg() ) );
    connect( ui->lineEdit_img, SIGNAL(textChanged(QString)),      this, SLOT( viewImg(QString) ) );
    connect( ui->sB_intensify, SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->pB_readGL,    SIGNAL( clicked() ),               this, SLOT( readGL() ));
    connect( ui->pB_saveBG,    SIGNAL( clicked() ),               this, SLOT( saveBG() ) );
    connect( ui->pB_saveOneSun,SIGNAL( clicked() ),               this, SLOT( saveOS() ) );
    connect( ui->pB_saveFlux,  SIGNAL( clicked() ),               this, SLOT( saveFX() ) );
    connect( ui->spinBox_i0,   SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->spinBox_j0,   SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->checkBox,     SIGNAL(stateChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->pB_ROI,       SIGNAL( clicked() ),               this, SLOT( setROI() ) );
    connect( ui->sB_iMin,      SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->sB_iMax,      SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->sB_jMin,      SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->sB_jMax,      SIGNAL(valueChanged(int)),         this, SLOT( setIntensify() ) );
    connect( ui->dSB_OS,       SIGNAL(valueChanged(double)),      this, SLOT( setIntensify() ) );
    connect( ui->pB_add,       SIGNAL( clicked() ),               this, SLOT( build() ) );

    namedWindow( "Flux Map", WINDOW_NORMAL );
    namedWindow( "Composite image", WINDOW_NORMAL );
    //Mat matDisplay(100,100,CV_64F);//to avoid warnings
    imshow( "Composite image",Cimg);
    imshow( "Flux Map", Cimg );
    setMouseCallback("Flux Map", mouseHandler, NULL);
}



ImageWorkBench::~ImageWorkBench()
{
    delete ui;
}



void ImageWorkBench::selectDir(){
    dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
                                                    "/run/media/marco/d2quadra/FluxMapper",
                                                    QFileDialog::ShowDirsOnly);
    ui->lineEdit_dir -> setText(dir);
}


void ImageWorkBench::selectImg(){
    dir=ui->lineEdit_dir -> text();
    filename = QFileDialog::getOpenFileName(this,tr("Open Image"),
                                                    dir,
                                                    tr("Image Files (*.png *.jpg *.bmp *.ppm)"));
    ui->lineEdit_img -> setText(filename);
}


void ImageWorkBench::setIntensify(){
    viewImg(filename);
}


void ImageWorkBench::viewImg(QString file){
    filename=file;
    img = imread(file.toStdString(),IMREAD_UNCHANGED);//CV_LOAD_IMAGE_UNCHANGED
    if(! img.data){
        cout << file.toStdString() <<" is not accessible!!!"<<"\n";
        return;
    }
    int height = img.rows;
    int width = img.cols;
    int chan = img.channels();
    int depth = img.depth();
    ui->lineEdit_width ->setText(QString::number(width));
    ui->lineEdit_height ->setText(QString::number(height));
    ui->comboBox_DataType ->setCurrentIndex(depth);
    ui->lineEdit_channel->setText(QString::number(chan));
    if(chan>1)
        cvtColor(img, img, COLOR_BGR2GRAY);
    int intensify = ui->sB_intensify -> value();
    Mat matDisplay(height,width,CV_16UC3);
    uint16_t value;
    int Value;
    double dval;
    double magn;
    if(depth==0 || depth==1)
        magn=pow(2.,16)/pow(2.,8);
    else if(depth==2 || depth==3)
        magn=pow(2.,16)/pow(2.,16);
    else if(depth==4 || depth==5)
        magn=pow(2.,16)/pow(2.,32);
    else
        magn=pow(2.,16)/pow(2.,64);
    for(int i=0; i<height; i++){
        for(int j=0; j< width; j++){
            if(depth==0 || depth==1)
                value=img.at<uint8_t>(i,j);
            else if(depth==2 || depth==3)
                value=img.at<uint16_t>(i,j);
            else if(depth==4 || depth==5)
                value=img.at<uint32_t>(i,j);
            else
                value=img.at<float_t>(i,j);
            Value=value*magn*intensify;
            if(Value < pow(2.,16)-1){
                matDisplay.at<Vec3s>(i,j)[0]=int(double(Value));
                matDisplay.at<Vec3s>(i,j)[1]=int(double(Value));
                matDisplay.at<Vec3s>(i,j)[2]=int(double(Value));
            }
            else{
                matDisplay.at<Vec3s>(i,j)[0]=0;
                matDisplay.at<Vec3s>(i,j)[1]=0;
                matDisplay.at<Vec3s>(i,j)[2]=short(65535);
            }
        }
    }
    matDisplay.copyTo(imgSt);
    cvtColor(imgSt, imgSt, COLOR_BGR2GRAY);
    Qt::CheckState state = ui->checkBox -> checkState();
    if(state==Qt::Checked){
        Point Pt,Pt1,Pt2;
        //draw screen border
        iMIN=ui->sB_iMin ->value();
        jMIN=ui->sB_jMin ->value();
        iMAX=ui->sB_iMax ->value();
        jMAX=ui->sB_jMax ->value();
        double Wscreen=ui->dSB_Wscreen->value();
        double Hscreen=ui->dSB_Hscreen->value();
        pxImDim=0.5*(Hscreen/double(iMAX-iMIN)+Wscreen/double(jMAX-jMIN));
        ui->dSB_pxImDim ->setValue(pxImDim);
        Pt1.x=jMIN;
        Pt1.y=iMIN;
        Pt2.x=jMAX;
        Pt2.y=iMAX;ui->sB_iMin ->setValue(iMIN);
        rectangle(matDisplay, Pt1, Pt2, Scalar(0,65535,0), 1);
        //draw anuli
        int j0=ui->spinBox_j0 -> value();
        int i0=ui->spinBox_i0 -> value();
        Pt.x=j0;
        Pt.y=i0;
        circle(matDisplay, Pt, 10, Scalar(0,65535,0), -1);
        int Nring=ui->spinBox_Nring -> value();
        int mmStep=ui->spinBox_step -> value();//step in mm
        ui->lineEdit_diameter->setText(QString::number(2*Nring*mmStep));
        double rStep=double(mmStep)/pxImDim; //step in px
        //cout << "rStep= "<<rStep<<"\n";
        int Radius=0;
        for(int i=1; i<=Nring;i++){
            Radius=int(double(i)*rStep+0.5);
            circle(matDisplay, Pt, Radius, Scalar(0,65535,0), 1);
        }
        //evaluate Sun*m^2
        int iR;
        double d;
        double glBG=ui->dSB_bg->value();
        double glOS=ui->dSB_OS->value();
        for(int i=0;i<100;i++){
            for(int j=0; j<4; j++)
                SumFlux[i][j]=0.;
        }
        for(int i=iMIN; i<iMAX; i++){
            for(int j=jMIN; j< jMAX; j++){
               value=img.at<uint16_t>(i,j);
               dval=Cimg.at<double>(i,j);
               d=sqrt(double(i-i0)*double(i-i0)+double(j-j0)*double(j-j0));
               iR=int(d/rStep);
               SumFlux[iR][0]=SumFlux[iR][0]+(double(value)-glBG)/(glOS);
               SumFlux[iR][2]=SumFlux[iR][2]+dval/glOS;//Cimg is background-free by definition
            }
        }
        QString fileFluxIntegral=dir+"/fluxIntegral.txt";
        QFile file2(fileFluxIntegral);
        if (!file2.open(QIODevice::WriteOnly | QIODevice::Text))
                return;
        QTextStream out(&file2);
        cout << "D(mm)\tFlux*Sanulus\tFluxTot(Suns*m^2)\tCimg_Flux*Sanulus\tCimg_FluxTot(Suns*m^2)"<<"\n";
        out <<  "D(mm)\tFlux*Sanulus\tFluxTot(Suns*m^2)\tCimg_Flux*Sanulus\tCimg_FluxTot(Suns*m^2)"<<"\n";
        for(int i=0;i<Nring;i++){
            SumFlux[i][0]=SumFlux[i][0]*pxImDim*pxImDim/1.E6;
            SumFlux[i][2]=SumFlux[i][2]*pxImDim*pxImDim/1.E6;
            if(i==0){
                SumFlux[i][1]=SumFlux[i][0];
                SumFlux[i][3]=SumFlux[i][2];
            }
            else{
                SumFlux[i][1]=SumFlux[i-1][1]+SumFlux[i][0];
                SumFlux[i][3]=SumFlux[i-1][3]+SumFlux[i][2];
            }
            cout << 2*mmStep*(i+1) <<"\t"<< SumFlux[i][0]<<"\t"<< SumFlux[i][1] <<"\t"<< SumFlux[i][2]<<"\t"<< SumFlux[i][3]<<"\n";
            out <<  2*mmStep*(i+1) <<"\t"<< SumFlux[i][0]<<"\t"<< SumFlux[i][1] <<"\t"<< SumFlux[i][2]<<"\t"<< SumFlux[i][3]<<"\n";
        }
        file2.close();
    }
    imshow("Flux Map",matDisplay);
}


void ImageWorkBench::readGL(){
    QString fileGrayLevel=dir+"/grayLevel.txt";
    QFile file3(fileGrayLevel);
    if (!file3.open(QIODevice::WriteOnly | QIODevice::Text))
            return;
    QTextStream out(&file3);
    int valMin=pow(2.,16);
    int valMax=0;
    double valMean=0.,valMeanQ=0;
    uint16_t value;
    double Npxs;
    cout << "Please select the ROI from Lup to Rdw"<<"\n";
    drag=0;
    while(drag!=2)
        waitKey(10);
    Npxs=double((mpxX2-mpxX1)*(mpxY2-mpxY1));
    printf("ROI (i %d ,j %d)Lup\t(i %d ,j %d)Rdw\tDi=%d\tDj=%d Npxs=%f\n",
           mpxY1,mpxX1,mpxY2,mpxX2,(mpxY2-mpxY1),(mpxX2-mpxX1),Npxs);
    out << mpxY1 <<"\t"<<mpxX1<<"\n";
    out << mpxY2 <<"\t"<<mpxX2<<"\n";
    for(int j=mpxX1; j< mpxX2; j++){
        for(int i=mpxY1; i<mpxY2; i++){
            value=img.at<uint16_t>(i,j);
            valMin=MIN(valMin,value);
            valMax=MAX(valMax,value);
            valMean=valMean+double(value)/Npxs;
            valMeanQ=valMeanQ+double(value)*double(value)/Npxs;
            if(i==mpxY1)
                out << j <<"\t";
            out << value;
            if(i<mpxY2-1)
                out<<"\t";
            else
                out<<"\n";
        }
    }
    file3.close();
    ui->lineEdit_mean -> setText(QString::number(int(valMean)));
    ui->lineEdit_sigma-> setText(QString::number(sqrt(valMeanQ-valMean*valMean)));
    ui->lineEdit_min  -> setText(QString::number(valMin));
    ui->lineEdit_max  -> setText(QString::number(valMax));
}

void ImageWorkBench::saveBG(){
    double gl=ui->lineEdit_mean -> text().toDouble();
    ui->dSB_bg -> setValue(gl);
    ui->dSB_AddBack-> setValue(gl);
}


void ImageWorkBench::saveOS(){
    double gl=ui->lineEdit_mean -> text().toDouble();
    double BG=ui->dSB_bg->value();
    gl=(gl-BG);
    ui->dSB_OS -> setValue(gl);
}

void ImageWorkBench::saveFX(){
    double gl=ui->lineEdit_max -> text().toDouble();
    double BG=ui->dSB_bg->value();
    double glOS=ui->dSB_OS->value();
    double gain=ui->dSB_addGain->value();
    gl=gl*gain;
    gl=(gl-BG);
    ui->dSB_Flux -> setValue(gl);
    gl=gl/glOS;
    ui->dSB_C -> setValue(gl);
}

void ImageWorkBench::build(){
    uint16_t value;
    int height=ui->lineEdit_height -> text().toInt();
    int width= ui->lineEdit_width  -> text().toInt();
    double BG=ui->dSB_bg ->value();
    double Again= ui->dSB_addGain ->value();
    BG=BG/Again;
    ui->dSB_AddBack->setValue(BG);
    double dval;
    for(int i=0; i<height; i++){
        for(int j=0; j< width; j++){
            value=img.at<uint16_t>(i,j);
            if(int(Again+0.5)==1){//new image!
                if(value <pow(2.,14)-1. )
                    if(double(value)-BG > 0.)
                        Cimg.at<double>(i,j)=double(value)-BG;
                    else
                        Cimg.at<double>(i,j)=0.;
                else
                    Cimg.at<double>(i,j)=-1.;
            }
            else{//add the new value if pixel_value is negative!
                dval=Cimg.at<double>(i,j);
                if(dval<0. && value <pow(2.,14)-1.){
                    if(double(value)-BG > 0.)
                        Cimg.at<double>(i,j)=(double(value)-BG)*Again;
                    else
                        Cimg.at<double>(i,j)=0.;
                }
            }
        }
    }
    imshow("Composite image",Cimg);
    waitKey(10);
}


void ImageWorkBench::setROI(){
    cout << "Please select the ROI from Lup to Rdw"<<"\n";
    drag=0;
    while(drag!=2)
        waitKey(10);
    printf("ROI (i %d ,j %d)Lup\t(i %d ,j %d)Rdw\tDi=%d\tDj=%d\n",
           mpxY1,mpxX1,mpxY2,mpxX2,(mpxY2-mpxY1),(mpxX2-mpxX1));
    iMIN=mpxY1;
    jMIN=mpxX1;
    iMAX=mpxY2;
    jMAX=mpxX2;
    ui->sB_iMin ->setValue(iMIN);
    ui->sB_jMin ->setValue(jMIN);
    ui->sB_iMax ->setValue(iMAX);
    ui->sB_jMax ->setValue(jMAX);
    double Wscreen=ui->dSB_Wscreen->value();
    double Hscreen=ui->dSB_Hscreen->value();
    pxImDim=0.5*(Hscreen/double(iMAX-iMIN)+Wscreen/double(jMAX-jMIN));
    ui->dSB_pxImDim ->setValue(pxImDim);
}

void ImageWorkBench::compute(){
    iMIN=ui->sB_iMin ->value();
    jMIN=ui->sB_jMin ->value();
    iMAX=ui->sB_iMax ->value();
    jMAX=ui->sB_jMax ->value();
    double Wscreen=ui->dSB_Wscreen->value();
    double Hscreen=ui->dSB_Hscreen->value();
    pxImDim=0.5*(Hscreen/double(iMAX-iMIN)+Wscreen/double(jMAX-jMIN));
    ui->dSB_pxImDim ->setValue(pxImDim);
    int j0=ui->spinBox_j0 -> value();
    int i0=ui->spinBox_i0 -> value();
    int Ncell=ui->sB_Ncell->value();
    int NcHalf=int(Ncell/2);
    double L=ui->dSB_Lmatrix->value();
    double cellSide=L/Ncell;
    double di,dj;
    double matrix[201][201][2]={0.};
    //uint16_t value;
    double dval;
    //double glBG=ui->dSB_bg->value();
    double glOS=ui->dSB_OS->value();
    int im,jm;
    cout << "Evaluating fluxMatrix ...";
    for(int i=iMIN; i<iMAX; i++){
        for(int j=jMIN; j< jMAX; j++){
           di=double(i-i0)*pxImDim;
           dj=double(j-j0)*pxImDim;
           if(fabs(di)<L/2. && fabs(dj)<L/2.){
               //value=img.at<uint16_t>(i,j);
               dval=Cimg.at<double>(i,j);
               im=NcHalf+int(di/cellSide);
               jm=NcHalf+int(dj/cellSide);
               if(0<=im && im<Ncell && 0<=jm && jm<Ncell){
                matrix[im][jm][0]++;
                matrix[im][jm][1]=matrix[im][jm][1]+dval/glOS;
               }
           }
        }
    }
    cout << "done!\n";
    QString fileMatrix=dir+"/fluxMatrix.txt";
    QFile file3(fileMatrix);
    if (!file3.open(QIODevice::WriteOnly | QIODevice::Text))
            return;
    QTextStream ffm(&file3);
    double flux;
    for(int i=0;i<Ncell;i++){
        for(int j=0;j<Ncell;j++){
            flux=0.;
            if(matrix[i][j][0]>0.)
                flux=matrix[i][j][1]/matrix[i][j][0];
            else
                flux=0.;
            if(j!=Ncell-1)
                ffm << flux <<"\t";
            else
                ffm << flux << "\n";
        }
    }
    file3.close();

    /*
    int nimg[3][2],height,width;
    nimg[0][0]=ui->sB_backMin -> value();
    nimg[0][1]=ui->sB_backMax -> value();
    nimg[1][0]=ui->sB_1sMin -> value();
    nimg[1][1]=ui->sB_1sMax -> value();
    nimg[2][0]=ui->sB_fluxMin -> value();
    nimg[2][1]=ui->sB_fluxMax -> value();
    QString fext,filename;
    Mat tmp;
    uint16_t value;
    //check image dimension
    filename=dir+"/img_Back_"+QString::number(nimg[0][0])+".ppm";
    tmp=imread(filename.toStdString(),CV_LOAD_IMAGE_UNCHANGED);
    if(! tmp.data){
        cout << filename.toStdString() <<" is not accessible!!!"<<"\n";
        return;
    }
    height = tmp.rows;
    width = tmp.cols;
    ui->lineEdit_width ->setText(QString::number(width));
    ui->lineEdit_height ->setText(QString::number(height));
    int arr[3] = {height,width,6};
    Mat L(3, arr, CV_64F, Scalar::all(0.0));
    double Ndata;
    //load and average images
    for(int jj=0;jj<3;jj++){
        if(jj==0)
            fext="Back";
        else if(jj==1)
            fext="OneS";
        else if(jj==2)
            fext="Flux";
        Ndata=double(nimg[jj][1]-nimg[jj][0]+1);
        for(int ii=nimg[jj][0];ii<=nimg[jj][1];ii++){
            filename=dir+"/img_"+fext+"_"+QString::number(ii)+".ppm";
            tmp=imread(filename.toStdString(),CV_LOAD_IMAGE_UNCHANGED);
            if(! tmp.data){
                cout << filename.toStdString() <<" is not accessible!!!"<<"\n";
                return;
            }
            for(int i=0;i<height;i++){
                for(int j=0;j<width;j++){
                    value=tmp.at<uint16_t>(i,j);
                    L.at<double>(i,j,jj*2)=L.at<double>(i,j,jj*2)+double(value)/Ndata;
                    L.at<double>(i,j,jj*2+1)=L.at<double>(i,j,jj*2+1)+double(value*value)/Ndata;
                }
            }
        }

    }
    double vDouble,vMean=0.,vMin=+99000.,vMax=-99000,vErr,vEmean=0;
    double B,eB,S,eS,F,eF;
    Mat map(height,width,CV_64F);
    Mat ERRmap(height,width,CV_64F);
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            B=L.at<double>(i,j,0);
            eB=L.at<double>(i,j,1)-B*B;
            if(eB>0.)
                eB=sqrt(eB);
            else{
                eB=0.;
                cout << "eB<0!" <<"\n";
            }
            S=L.at<double>(i,j,2);
            eS=L.at<double>(i,j,3)-S*S;
            if(eS>0.)
                eS=sqrt(eS);
            else{
                eS=0.;
                cout << "eS<0!" <<"\n";
            }
            F=L.at<double>(i,j,4);
            eF=L.at<double>(i,j,5)-F*F;
            if(eF>0.)
                eF=sqrt(eF);
            else{
                eF=0.;
                cout << "eF<0!" <<"\n";
            }
            vDouble=(F-B)/(S-B);
            vErr=(eF+eB)/(F-B)+(eS+eB)/(S-B);
            map.at<double>(i,j)=vDouble;
            ERRmap.at<double>(i,j)=vDouble*vErr;
            vMean=vMean+vDouble;
            vMin=min(vMin,vDouble);
            vMax=max(vMax,vDouble);
            vEmean=vEmean+vDouble*vErr;
        }
    }
    imshow( "Flux Map", map );
    vMean=vMean/double(height*width);
    vEmean=vEmean/double(height*width);
    ui->lineEdit_mean ->setText(QString::number(vMean)+QChar(0xb1)+QString::number(vEmean));
    ui->lineEdit_min ->setText(QString::number(vMin));
    ui->lineEdit_max ->setText(QString::number(vMax));
*/
}



void mouseHandler(int event, int x, int y, int flags, void* param){
    /* user press left button */
    Point point;
    if (event == EVENT_LBUTTONDOWN){
        point = Point(x, y);
        mpxX1=x;
        mpxY1=y;
        drag=1;
        int value=imgSt.at<uint16_t>(x,y);
        printf("left button DW i=%d j=%d value=%d = %f\n",mpxY1,mpxX1,value,value/(pow(2.,16)-1));

    }
    if (event == EVENT_LBUTTONUP){
        point = Point(x, y);
        mpxX2=x;
        mpxY2=y;
        if(mpxX2 > mpxX1 && mpxY2 > mpxY1)
         drag=2;
        else
         drag=1;
        //printf("left button UP i=%d j=%d\n",mpxY2,mpxX2);
    }
    if (event == EVENT_RBUTTONDOWN) {
      drag=-1;
      //printf("right button\n");
    }
    if (event == EVENT_MBUTTONDOWN) {
      drag=10;
      //printf("click medium button\n");
    }

}
