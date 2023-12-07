/*Author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it

kVISdishis the master software for infield experimental evaluation of 3D-shape of
solar dishes according to the VIS approach: a suitably wide LCD monitor is
placed orthogonally and centrally to the paraboloid axis at the
point (optically) conjugated with the observation point, always placed along the paraboloid axis.
The LCD monitor is connected to the “slave” PC,
while the GigE camera, placed at the observation point, is connected to the neighbor “master” PC.
KVISdish and kVISdishSlave have to be installed on the “master” and “slave” PC, respectively.
“Master” and “slave” PCs are connected by Ethernet


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
#include "kVISdish.h"
#include "ui_kVISdish.h"
#include <QApplication>
#include <QtGui>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <QFileDialog>

using namespace std;
using namespace cv;

/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia
Roma - Italy*/


//functions
void unitV12(double x1, double y1, double z1, double x2, double y2, double z2);
void unitVnIde(double x, double y, double focal);
void pxColor(double val);
inline QImage  cvMatToQImage( const cv::Mat &inMat );
inline QPixmap cvMatToQPixmap( const cv::Mat &inMat );

//known everywhere
QString Path="/run/media/marco/d2quadra/VISdish/";
QString PathSlave="smb://192.168.138.73/marco/Immagini/blobPattern";
QString PathSyncGrab="/opt/Vimba_1_3/VimbaCPP/Examples/SynchronousGrab/Console/Build/Make/binary/x86_64bit/SynchronousGrabConsole /f:/run/media/marco/d2quadra/VISdish/frames/img";
int Nfacet=1,NblobPattern=0;
double xyf[100][2];
double Vservice[3];//unit vector
double caronte;



Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);

    //loading facet coordinates
    int svalue;
    Nfacet=0;
    QFile ff(Path+"XYfacet.txt");
    if(!ff.open(QIODevice::ReadOnly | QIODevice::Text)){
        printf("error opening XYfacet.txt!\n");
        return;
    }
    QTextStream stream(&ff);
    do{
        stream >> svalue >> xyf[Nfacet][0] >> xyf[Nfacet][1];
        //printf("%d %f %f \n",svalue,xyf[Nfacet][0],xyf[Nfacet][1]);
        Nfacet++;
    } while (!stream.atEnd());
    ff.close();
    Nfacet--;
    //printf("Nfacet = %d \n",Nfacet);
    ui -> sB_Nfacet -> setMaximum(Nfacet);
    ui -> sB_Nfacet -> setRange(1,Nfacet);
    ui -> sB_Nfacet -> setValue(Nfacet);
    ui -> dSB_Xmirror -> setValue(xyf[0][0]);
    ui -> dSB_Ymirror -> setValue(xyf[0][1]);
    ui -> lineEdit_slaveFolder -> setText(PathSlave);

    //signals&slots
    connect(ui-> pushButton_plotBlob,SIGNAL(clicked(bool)),this,SLOT(plotBlob()));
    connect(ui-> pushButton_TakePicture,SIGNAL(clicked(bool)),this,SLOT(takeOnePicture()));
    connect(ui-> pushButton_process,SIGNAL(clicked(bool)),this,SLOT(Process()));
}



Widget::~Widget()
{
    delete ui;
}



void Widget::on_pushButton_slaveFolder_clicked()
{
    PathSlave=QFileDialog::getExistingDirectory(
                this,
                "Choose the directory for VISdish-slave where save the blob pattern", //titolo della finestra
                Path, //directory iniziale
                QFileDialog::ShowDirsOnly); //tipi di file da cercare
    ui-> lineEdit_slaveFolder -> setText(PathSlave);
}



void Widget::plotBlob()
{
    int width,hight,rBlob,iBW,iOk,xBlob,yBlob,kmin,kmax,iB1L2;
    double Lwidth,Lhight;
    Point2f Pt1,Pt2,Pt3,Pt4;
    Mat mat;
    double focal,zObs,zMonitor,x,y,z,vi[3],vr[3],vn[3],cosinc,t,xs,ys;
    QString Line,filename;
    Qt::CheckState state;
    cv::Scalar scalar;

    //check if PathSlave was set
    iOk=1;
    Line= ui-> lineEdit_slaveFolder -> text();
    if(Line=="!!!!! please select the remote folder for VISdish-slave !!!!")
        iOk=0;
    else{
        PathSlave= ui-> lineEdit_slaveFolder -> text();
        printf("PathSlave = %s\n",(PathSlave.toStdString()).c_str());
    }

    //check&reset checkBox
    state = ui->checkBox -> checkState();
    if(state==Qt::Unchecked) ui -> checkBox -> setCheckState ( Qt::Checked );

    //read parameters
    width  = ui -> sB_MpxW  -> value();
    hight  = ui -> sB_MpxH  -> value();
    Lwidth = ui -> dSB_MmmW -> value();
    Lhight = ui -> dSB_MmmH -> value();
    rBlob  = ui -> sB_rBlob -> value();

    //Blob or Line?
    iB1L2 = ui-> comboBox_BlobLine -> currentIndex();
    iB1L2++;

    //color & mat
    iBW = ui -> comboBox -> currentIndex();
    if(iBW==0){
        scalar.val[0]=255;
        scalar.val[1]=255;
        scalar.val[2]=255;
        mat = Mat::zeros(hight, width, CV_8UC1);
    }
    else{
        scalar.val[0]=0;
        scalar.val[1]=0;
        scalar.val[2]=0;
        mat = Mat::ones(hight, width, CV_8UC1)*255;
    }

    //Plot All Points
    state = ui->checkBox_PAP-> checkState();
    if(state==Qt::Checked) {
     kmin= ui -> sB_Nfacet -> minimum();
     kmax= ui -> sB_Nfacet -> maximum();
    }
    else{
     Nfacet = ui -> sB_Nfacet -> value();
     kmin=Nfacet;
     kmax=Nfacet;
    }
    kmin--;
    kmax--;
    focal= ui -> dSB_focal -> value();
    zObs= ui -> dSB_observer -> value();
    zMonitor = ui -> dSB_zMonitor -> value();

    QFile fileCTRL(PathSlave+"/currentBlob.ctrl");
    fileCTRL.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&fileCTRL);
    state = ui->checkBox_FS -> checkState();
    if(state==Qt::Unchecked)
        out << "normalscreen";
    else
        out << "fullscreen";
    out << "\n";
    out << hight <<"\t"<< width <<"\n";

    for(int k=kmin;k<=kmax;k++){
        x=xyf[k][0];
        y=xyf[k][1];
        unitVnIde(x,y,focal);
        for(int i=0;i<3;i++)
            vn[i]=Vservice[i];
        z=caronte;
        unitV12(0., 0., zObs, x, y, z);
        cosinc=0.;
        for(int i=0;i<3;i++){
            vi[i]=Vservice[i];
            cosinc=cosinc+vi[i]*vn[i];
        }
        cosinc=-cosinc;
        for(int i=0;i<3;i++){
           vr[i]=vi[i]+2.0*cosinc*vn[i];
        }
        t=(zMonitor-z)/vr[2];
        xs=(x+vr[0]*t)*1000;//mm
        ys=(y+vr[1]*t)*1000;//mm
        ui -> dSB_Xmirror -> setValue(xyf[k][0]);
        ui -> dSB_Ymirror -> setValue(xyf[k][1]);
        ui -> dSB_Xblob ->setValue(xs);
        ui -> dSB_Yblob ->setValue(ys);
//    printf("vi: %f %f %f\n",vi[0],vi[1],vi[2]);
//    printf("vn: %f %f %f\n",vn[0],vn[1],vn[2]);
//    printf("vr: %f %f %f\n",vr[0],vr[1],vr[2]);
//    printf("zObs = %f alpha = %f t = %f cosinc = %f\n",zObs,alpha,t,cosinc);
        xBlob=int(xs/(Lwidth/2.)*(width/2)+width/2);
        yBlob=int(ys/(Lhight/2.)*(hight/2)+hight/2);
        if(iB1L2==1){
            Pt1.x= xBlob;
            Pt1.y= yBlob;
            circle(mat, Pt1, rBlob, scalar, -1);
            out << iB1L2 <<"\t"<< iBW <<"\t"<< rBlob <<"\n";
            out << Pt1.x <<"\t"<< Pt1.y <<"\n";
        }
        else if(iB1L2==2){
            Pt1.x= xBlob;
            Pt1.y= 0;
            Pt2.x= xBlob;
            Pt2.y= hight;
            Pt3.x= 0;
            Pt3.y= yBlob;
            Pt4.x= width;
            Pt4.y= yBlob;
            line(mat,Pt1,Pt2,scalar,2*rBlob,8,0);
            line(mat,Pt3,Pt4,scalar,2*rBlob,8,0);
            out << iB1L2 <<"\t"<< iBW  <<"\t"<< rBlob <<"\n";
            out << Pt1.x <<"\t"<< Pt1.y <<"\t"<< Pt2.x <<"\t"<< Pt2.y <<"\n";
            out << iB1L2 <<"\t"<< iBW  <<"\t"<< rBlob <<"\n";
            out << Pt3.x <<"\t"<< Pt3.y <<"\t"<< Pt4.x <<"\t"<< Pt4.y <<"\n";
        }
    }

    imwrite(Path.toStdString()+"currentBlob.ppm",mat);
    if(iOk==1){
        filename = PathSlave+"/currentBlob.ppm";
        imwrite(filename.toStdString(),mat);
    }
    ui->label_10->setPixmap(cvMatToQPixmap(mat));
    fileCTRL.close();
    waitKey(10);//to ensure GUI reactivity
    mat.release();

}



void Widget::on_pushButton_clicked(){
    int width,hight,rBlob,nBlob,npx,iBW,Dj,Di,iB1L2,err;
    double Lwidth,Lhight,step;
    Point2f Pt1,Pt2;
    Mat mat;
    QString Line,filename,command;
    Qt::CheckState state;
    cv::Scalar scalar;

    //check if PathSlave was set
    Line= ui-> lineEdit_slaveFolder -> text();
    if(Line=="!!!!! please select the remote folder for VISdish-slave !!!!"){
        printf("!!!!! please select the remote folder for VISdish-slave !!!!\n");
        return;
    }

    //check&reset checkBox
    state = ui->checkBox -> checkState();
    if(state==Qt::Unchecked) ui -> checkBox -> setCheckState ( Qt::Checked );

    //read parameters
    width  = ui -> sB_MpxW  -> value();
    hight  = ui -> sB_MpxH  -> value();
    Lwidth = ui -> dSB_MmmW -> value();
    Lhight = ui -> dSB_MmmH -> value();
    rBlob  = ui -> sB_rBlob -> value();
    nBlob  = ui -> sB_Nblob -> value();

    //Blob or Line?
    iB1L2 = ui-> comboBox_BlobLine -> currentIndex();
    iB1L2++;

    //color
    iBW = ui -> comboBox -> currentIndex();
    if(iBW==0){
        scalar.val[0]=255;
        scalar.val[1]=255;
        scalar.val[2]=255;
        mat = Mat::zeros(hight, width, CV_8UC1);
    }
    else{
        scalar.val[0]=0;
        scalar.val[1]=0;
        scalar.val[2]=0;
        mat = Mat::ones(hight, width, CV_8UC1)*255;
    }

    //
    npx=hight;
    Dj=int((width-hight)/2);
    Di=0;
    if(width<hight) {
        npx=width;
        Di=int((hight-width)/2);
        Dj=0;
    }
    step=double(npx-2*rBlob)/double(nBlob-1);

    //open and write header of the log file
    QFile file(Path+"infoBlob.txt");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
         printf("file non accessibile!\n");
         return;
     }
    QTextStream stream ( &file );
    stream << "Monitor(mm):\t"<<Lhight<<"x"<<Lwidth<<"\n";
    stream << "Monitor(px):\t"<<hight<<"x"<<width<<"\n";
    int jMin,jMax;
    if(iB1L2 == 1){
        stream << "blobRadius:\t"<<rBlob<<"\n";
        stream << "blobNumber:\t"<<nBlob<<"x"<<nBlob<<"\n";
        stream << "kind:\tBlob"<<"\n";
    }
    else if(iB1L2==2){
        stream << "lineTichness:\t"<<2*rBlob<<"\n";
        stream << "lineNumber:\t"<<2<<"x"<<nBlob<<"\n";
        stream << "kind:\tLine"<<"\n";
    }
    //start loop

    NblobPattern=0;
    for(int i=0;i<nBlob;i++){
        if(iB1L2 == 1){
            jMin=0;
            jMax=nBlob;}
        else if(iB1L2==2){
            jMin=i;
            jMax=i+1;
        }
        for(int j=jMin;j<jMax;j++){
            if(iB1L2 == 1){ //blob
//blob
                NblobPattern++;
                //color & mat
                iBW = ui -> comboBox -> currentIndex();
                if(iBW==0){
                    scalar.val[0]=255;
                    scalar.val[1]=255;
                    scalar.val[2]=255;
                    mat = Mat::zeros(hight, width, CV_8UC1);
                }
                else{
                    scalar.val[0]=0;
                    scalar.val[1]=0;
                    scalar.val[2]=0;
                    mat = Mat::ones(hight, width, CV_8UC1)*255;
                }
                Pt1.x= int(rBlob+step*j)+Dj;
                Pt1.y= int(rBlob+step*i)+Di;
                QFile fileCTRL(PathSlave+"/currentBlob.ctrl");
                fileCTRL.open(QIODevice::WriteOnly | QIODevice::Text);
                QTextStream out(&fileCTRL);
                state = ui->checkBox_FS -> checkState();
                if(state==Qt::Unchecked)
                    out << "normalscreen";
                else
                    out << "fullscreen";
                out << "\n";
                out << hight <<"\t"<< width <<"\n";
                out << iB1L2 <<"\t"<< iBW <<"\t"<< rBlob <<"\n";
                out << Pt1.x <<"\t"<< Pt1.y  <<"\n";
                fileCTRL.close();
                circle(mat, Pt1, rBlob, scalar, -1);
                stream <<i<<"\t"<<j<<"\t"<<Pt1.x<<"\t"<<Pt1.y<<"\t"<<Pt1.x/width*Lwidth<<"\t"<<Pt1.y/hight*Lhight<<"\n";
                state = ui-> checkBox_save -> checkState();
                if(state==Qt::Checked) {
                    filename = Path +  "blob_" + QString::number(NblobPattern) +".ppm";
                    imwrite(filename.toStdString(),mat);
                }
                ui->label_10->setPixmap(cvMatToQPixmap(mat));
                waitKey(10);//to ensure GUI reactivity
                command=PathSyncGrab+"_"+QString::number(NblobPattern);
                err=system((command.toStdString()).c_str());
                if(err!=0) cout<<"system returned an error!!!";
            }
            else{ //line
//horizontal line
                NblobPattern++;//horizontal line
                //color & mat
                iBW = ui -> comboBox -> currentIndex();
                if(iBW==0){
                    scalar.val[0]=255;
                    scalar.val[1]=255;
                    scalar.val[2]=255;
                    mat = Mat::zeros(hight, width, CV_8UC1);
                }
                else{
                    scalar.val[0]=0;
                    scalar.val[1]=0;
                    scalar.val[2]=0;
                    mat = Mat::ones(hight, width, CV_8UC1)*255;
                }
                Pt1.x= Dj;
                Pt1.y= int(rBlob+step*i)+Di;
                Pt2.x= width-Dj;
                Pt2.y= int(rBlob+step*i)+Di;
                QFile fileCTRL1(PathSlave+"/currentBlob.ctrl");
                fileCTRL1.open(QIODevice::WriteOnly | QIODevice::Text);
                QTextStream out1(&fileCTRL1);
                state = ui->checkBox_FS -> checkState();
                if(state==Qt::Unchecked)
                    out1 << "normalscreen";
                else
                    out1 << "fullscreen";
                out1 << "\n";
                out1 << hight <<"\t"<< width <<"\n";
                out1 << iB1L2 <<"\t"<< iBW <<"\t"<< rBlob <<"\n";
                out1 << Pt1.x <<"\t"<< Pt1.y <<"\t"<< Pt2.x <<"\t"<< Pt2.y <<"\n";
                fileCTRL1.close();
                line(mat,Pt1,Pt2,scalar,2*rBlob,8,0);
                stream <<i<<"\t"<<"hori"<<"\t"<<
                        Pt1.x<<"\t"<<Pt1.y<<"\t"<<Pt1.x/width*Lwidth<<"\t"<<Pt1.y/hight*Lhight<<"\t"<<
                        Pt2.x<<"\t"<<Pt2.y<<"\t"<<Pt2.x/width*Lwidth<<"\t"<<Pt2.y/hight*Lhight<<"\n";
                state = ui-> checkBox_save -> checkState();
                if(state==Qt::Checked) {
                    filename = Path +  "line_" + QString::number(NblobPattern) +".ppm";
                    imwrite(filename.toStdString(),mat);
                }
                ui->label_10->setPixmap(cvMatToQPixmap(mat));
                waitKey(10);//to ensure GUI reactivity
                command=PathSyncGrab+"_"+QString::number(NblobPattern);
                err=system((command.toStdString()).c_str());
                if(err!=0) cout<<"system returned an error!!!";
//vertical line
                NblobPattern++;//vertical line
                //color & mat
                iBW = ui -> comboBox -> currentIndex();
                if(iBW==0){
                    scalar.val[0]=255;
                    scalar.val[1]=255;
                    scalar.val[2]=255;
                    mat = Mat::zeros(hight, width, CV_8UC1);
                }
                else{
                    scalar.val[0]=0;
                    scalar.val[1]=0;
                    scalar.val[2]=0;
                    mat = Mat::ones(hight, width, CV_8UC1)*255;
                }
                Pt1.x= int(rBlob+step*j)+Dj;
                Pt1.y= 0;
                Pt2.x= int(rBlob+step*j)+Dj;
                Pt2.y= hight;
                QFile fileCTRL2(PathSlave+"/currentBlob.ctrl");
                fileCTRL2.open(QIODevice::WriteOnly | QIODevice::Text);
                QTextStream out2(&fileCTRL2);
                state = ui->checkBox_FS -> checkState();
                if(state==Qt::Unchecked)
                    out2 << "normalscreen";
                else
                    out2 << "fullscreen";
                out2 << "\n";
                out2 << hight <<"\t"<< width <<"\n";
                out2 << iB1L2 <<"\t"<< iBW <<"\t"<< rBlob <<"\n";
                out2 << Pt1.x <<"\t"<< Pt1.y <<"\t"<< Pt2.x <<"\t"<< Pt2.y <<"\n";
                fileCTRL2.close();
                line(mat,Pt1,Pt2,scalar,2*rBlob,8,0);
                stream <<i<<"\t"<<"vert"<<"\t"<<
                        Pt1.x<<"\t"<<Pt1.y<<"\t"<<Pt1.x/width*Lwidth<<"\t"<<Pt1.y/hight*Lhight<<"\t"<<
                        Pt2.x<<"\t"<<Pt2.y<<"\t"<<Pt2.x/width*Lwidth<<"\t"<<Pt2.y/hight*Lhight<<"\n";
                state = ui-> checkBox_save -> checkState();
                if(state==Qt::Checked) {
                    filename = Path +  "line_" + QString::number(NblobPattern) +".ppm";
                    imwrite(filename.toStdString(),mat);
                }
                ui->label_10->setPixmap(cvMatToQPixmap(mat));
                waitKey(10);//to ensure GUI reactivity
                command=PathSyncGrab+"_"+QString::number(NblobPattern);
                err=system((command.toStdString()).c_str());
                if(err!=0) cout<<"system returned an error!!!";
            }
            state = ui->checkBox -> checkState();
            if(state==Qt::Unchecked) {
               i=nBlob;
               j=jMax;
            }
        }
    }
    file.close();
//    destroyAllWindows();
    mat.release();

}


void Widget::takeOnePicture(){
    QString command=PathSyncGrab+"_shot";
    int err=system((command.toStdString()).c_str());
    if(err!=0) cout<<"system returned an error!!!";
}


void Widget::Process(){
    QString PathPro="/run/media/marco/d2quadra/VISdish/misura6nov2015/frames/";
    Mat matOriginal,matCropped,mat,dst;
    //open final result file
    QFile file(PathPro+"results4.txt");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)){
         printf("file non accessibile!\n");
         return;
    }
    QTextStream stream ( &file );
    int Nframe=270;
    int thres1= ui -> spinBox_threshold -> value();   //threshold
    int thres0=20;//used to delimit the dish boundaries
    int jROI=240;   //ROI
    int iROI=30;    //ROI
    int wROI=1960;  //ROI
    int hROI=1960;  //ROI
    int reduction=4;//reduction of resolution
    double factorRed= 1./double(reduction);
    int width=wROI/reduction;  //final width
    int height=hROI/reduction; //final height
    double focal=7016; //dish focal (mm) experimental value by FluxMapper in Moon tracking
    double dishDiam=11554.; //dish diameter (mm) misura di Pino del 16 Novembre 2015
    double d=319000;   //observer distance (mm)
    double dRim=d-0.25/focal*0.25*dishDiam*dishDiam;//rim-observer distance
    double fpx=double(width)*dRim/dishDiam;//camera focal length in pixel obtained by dishDiam/dRim=width/fpx
    double tiltX=0.01242; //correction for aiming error
    double tiltY=0.00017;//idem
    namedWindow("image", WINDOW_NORMAL);
 //   mat=imread(PathPro.toStdString()+"img_0",CV_LOAD_IMAGE_ANYDEPTH);
 //   int width    = mat.cols;
 //   int height   = mat.rows;
 //   int depth    = mat.depth();
    printf("Image %dx%d\n",width,height);
    int arr[3] = {width,height,5};
    Mat L(3, arr, CV_64FC1, Scalar::all(0));
    //Mat L(3, arr, CV_16UC(1), Scalar::all(0));
    //L(i,j,0)  if( img0(i,j) is under treshold)
    //             0
    //          else
    //             N frames over threshold; in the end is used for storaging slopeDev
    //L(i,j,1)  NframeMIN for horizontal line
    //L(i,j,2)  NframeMAX for horizontal line
    //L(i,j,3)  NframeMIN for vertical line
    //L(i,j,4)  NframeMAX for vertical line
    Mat sampled=Mat::zeros(width,height,CV_8UC1);
    QString fileImg;
    int intensity,h0v1,iM,thres;
    printf("\n>>>>> processing is started ...\n");
    fflush(stdout);
    for(int n=0;n<=Nframe;n++){
        thres=thres1;
        if(n==0) thres=thres0;
        if(2*int(double(n)/2.)-n < 0)
            h0v1=0; //horizontal line
        else
            h0v1=1; //vertical line
//        printf("n=%d h1v2=%d\n",n,h1v2);
        fileImg=PathPro+"img_"+QString::number(n);//image to be loaded
        matOriginal=imread(fileImg.toStdString(),IMREAD_ANYDEPTH);//original image
        matCropped = matOriginal(Rect(jROI,iROI,wROI,hROI)); //crop image with ROI
        cv::resize(matCropped, mat, Size(), factorRed, factorRed, INTER_AREA);//reduce image resolution of factorRed
//        imwrite("/home/marco/Immagini/matResized.jpg",mat);
        threshold(mat,dst,thres,254,THRESH_BINARY ); //thresholded image
        imshow("image",dst);
        waitKey(10);
        if(n==0) waitKey(2000);
        for(int i=0;i<height;i++){
            for(int j=0;j<width;j++){
                intensity= mat.at<uchar>(i, j);
//                printf("%d %d %d: intensity = %d\n",n,j,i,intensity);
                if(intensity>=thres){
                    if(n==0)
                        L.at<double>(i,j,0)=1.;
                    else{
                        if(L.at<double>(i,j,0)>0.){
                            L.at<double>(i,j,0)++;
                            iM=1;
                            if(L.at<double>(i,j,iM+2*h0v1)>0.) iM=2;
                            L.at<double>(i,j,iM+2*h0v1)=double(n);
                     }
                    }
                }
            }
        }
//        printf("Nframe = %d\n",n);
    }
    double count=0.;
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            if(L.at<double>(i,j,0)>1. && L.at<double>(i,j,1) > 0. && L.at<double>(i,j,3) > 0.){
                count++;
                sampled.at<uchar>(i,j)=255;
//                printf("L_%d_%d: %d %d %d %d %d\n",j,i,L.at<double>(i,j,0),L.at<double>(i,j,1),L.at<double>(i,j,2),L.at<double>(i,j,3),L.at<double>(i,j,4));

            }
        }
    }
    namedWindow("sampledSurface", WINDOW_NORMAL);
    imshow("sampledSurface",sampled);
    imwrite((PathPro+"sampledSurface.jpg").toStdString(),sampled);
    waitKey(10);
    double fracSamArea=count/(double(height*width)*3.14/4.);
    printf("image processing completed!\n The fraction of sampled area is about %f\n",fracSamArea);
    fflush(stdout);

    int Nh,Nh2,Nv,Nv2,Dn,jSqua,iSqua;
    double xj,yi,vr,vz,t,x,y,z,xsr,ysr,xs,ys,zs,modulo,costeta,
           slopeDev,sum,sum2,slopeDevX,sumX,sumX2,slopeDevY,sumY,sumY2,
           sumDn,sdMAX;
    double uvi[3],uvn[3],uvr[3],uvnid[3];
    double zMon=7242.;// z monitor
    double Dxs=0.;//x-offset of monitor center (mm)
    double Dys=5.;//y-offset of monitor center (mm)
    double alpha=0.0034;//monitor rotation
    int wPx=1920;//width monitor (px)
    int hPx=1080;//height monitor (px)
    int Dj=(wPx-hPx)/2;//exceeding lateral band of the monitor
    int hhPx=hPx/2;//half height (px)
    double wMon=1020.;//width monitor (mm)
//    double hMon=570.;//height monitor (mm)
    double px2mm=wMon/wPx;
    int linThick=8;//line thickness (px)
    count=0.;
    sum=.0;
    sum2=.0;
    sumX=.0;
    sumX2=.0;
    sumY=.0;
    sumY2=.0;
    sumDn=0.;
    sdMAX=0.;
    sumX=0.;
    sumY=0.;
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            if(L.at<double>(i,j,1)>0. && L.at<double>(i,j,3)>0.){
                //coordinates of the point on the dish mirror surface
                count++;
                xj=+(double(j)-double(width)/2.); //the sign is inverted to straighten the immage
                yi=-(double(i)-double(height)/2.);//idem
                unitV12(xj,yi,fpx,0.,0.,0.);//reflected unit vector calculated in the camera pin-hole model)
                for(int ii=0;ii<3;ii++)
                    uvr[ii]=Vservice[ii];
                vr=sqrt(uvr[0]*uvr[0]+uvr[1]*uvr[1]);
                vz=uvr[2];
                t=(vz+sqrt(vz*vz+vr*vr*d/focal))/vr/vr*2.*focal;
                x=uvr[0]*t;  //coordinates of mirrorPoint
                y=uvr[1]*t;
                z=uvr[2]*t+d;
                unitVnIde(x,y,focal);//unit vector normal to the ideal paraboloid
                for(int ii=0;ii<3;ii++)
                    uvnid[ii]=Vservice[ii];
                // coordinate of the point source on monitor
                //Nframe of the involved horizontal line
                Nh=int(L.at<double>(i,j,1)+0.5);
                Dn=0;
                if(L.at<double>(i,j,2)>0.){
                    Nh2=int(L.at<double>(i,j,2)+0.5);
                    Dn=Nh2-Nh;
                    if(Dn>2){
                        Nh=(Nh+Nh2)/2;
                        if(2*int(double(Nh)/2.)-Nh == 0) Nh--;
                    }
                }
                sumDn=sumDn+Dn;
                //Nframe of the involved vertical line
                Nv=int(L.at<double>(i,j,3)+0.5);
                Dn=0;
                if(L.at<double>(i,j,4)>0.){
                    Nv2=int(L.at<double>(i,j,4)+0.5);
                    Dn=Nv2-Nv;
                    if(Dn>2){
                        Nv=(Nv+Nv2)/2;
                        if(2*int(double(Nv)/2.)-Nv < 0) Nv--;
                    }
                }
                sumDn=sumDn+Dn;
                //monitor coordinates
                jSqua=linThick/2+int(Nv/2)*linThick+Dj;//Nv/2 because nframe includes horizontal and vertical patterns
                iSqua=linThick/2+int(Nh/2)*linThick;   //Nh/2 idem
                xsr=jSqua-hhPx;//monitor is facing the dish!!!
                ysr=iSqua-hhPx;
                //square-source coordinate in dish frame (mm)
                xs=(xsr*px2mm+Dxs)*cos(alpha) +(ysr*px2mm+Dys)*sin(alpha);
                ys=-(xsr*px2mm+Dxs)*sin(alpha) +(ysr*px2mm+Dys)*cos(alpha);
                zs=zMon;
                //incidence unit-vector
                unitV12(xs,ys,zs,x,y,z);
                modulo=0.;
                for(int ii=0;ii<3;ii++){
                    uvi[ii]=Vservice[ii];
                    uvn[ii]=uvr[ii]+uvi[ii];// Vnormal is the vectorial sum Vinc + Vref
                    modulo=modulo+uvn[ii]*uvn[ii];
                }
                for(int ii=0;ii<3;ii++){
                    uvn[ii]=-uvn[ii]/sqrt(modulo);//- to point towards the observer
                }
                modulo=0.;
                uvn[0]=cos(acos(uvn[0])+tiltX);//mis-aiming correction
                uvn[1]=cos(acos(uvn[1])+tiltY);
                modulo=uvn[0]*uvn[0]+uvn[1]*uvn[1]+uvn[2]*uvn[2];
                for(int ii=0;ii<3;ii++)
                    uvn[ii]=uvn[ii]/sqrt(modulo);
                //computing slope deviation
                costeta=.0;
                for(int ii=0;ii<3;ii++)
                    costeta=costeta+uvn[ii]*uvnid[ii];
                slopeDev=acos(costeta)*1000.;//slope deviation (mrad)
                slopeDevX=(acos(uvn[0])-acos(uvnid[0]))*1000.;//slopeX deviation (mrad)
                slopeDevY=(acos(uvn[1])-acos(uvnid[1]))*1000.;//slopeY deviation (mrad)
                L.at<double>(i,j,0)=slopeDev;//store slopDev value
                sum=sum+slopeDev;
                sum2=sum2+slopeDev*slopeDev;
                sumX=sumX+slopeDevX;
                sumX2=sumX2+slopeDevX*slopeDevX;
                sumY=sumY+slopeDevY;
                sumY2=sumY2+slopeDevY*slopeDevY;
                sdMAX=max(sdMAX,slopeDev);
//                if(count>149736. && count <149750.){
//                    printf("Ndato= %f\n",count);
//                    printf("P(j=%d, i=%d): %f %f %f\n",j,i,x,y,z);
//                    printf("Nv=%d -> jSqua=%d Nh=%d -> iSqua=%d\n",Nv,jSqua,Nh,iSqua);
//                    printf("S: %f %f %f\n",xs,ys,zs);
//                    printf("uvn_exper: %f %f %f\n",uvn[0],uvn[1],uvn[2]);
//                    printf("uvn_ideal: %f %f %f\n",uvnid[0],uvnid[1],uvnid[2]);
//                    printf("slopDev (mrad)= %f\n",slopeDev);
//                    //fflush(stdout);
//                    //waitKey(0);
//                }
                stream << x <<"\t"<< y <<"\t"<< z <<"\t"<< uvn[0] <<"\t"<< uvn[1] <<"\t"<< uvn[2] <<"\t"<< slopeDev <<"\t"<< slopeDevX <<"\t"<< slopeDevY <<"\n";
            }
        }
    }
    // creation and plot of slope Deviation Map
    Mat slopeDevMap=Mat::zeros(width,height,CV_8UC3);
    double val;
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            if(L.at<double>(i,j,1)>0. && L.at<double>(i,j,3)>0.){
                slopeDev=L.at<double>(i,j,0);
                val=slopeDev/sdMAX;
                pxColor(val);
                slopeDevMap.at<Vec3b>(i,j)[0] = Vservice[0];
                slopeDevMap.at<Vec3b>(i,j)[1] = Vservice[1];
                slopeDevMap.at<Vec3b>(i,j)[2] = Vservice[2];
            }
        }
    }
    namedWindow("slopeDevMap", WINDOW_NORMAL);
    imshow("slopeDevMap",slopeDevMap);
    imwrite((PathPro+"slopeDevMap.jpg").toStdString(),slopeDevMap);
    double sdm=sum/count;
    double sigma_sdm=(sum2-sum*sum/count)/double(count-1);
    sigma_sdm=sqrt(sigma_sdm);
    double Dframe=double(sumDn)/count;
    printf("%d data have been succesfully processed!\n",int(count+0.5));
    printf("<DeltaFrame> = %f\n",Dframe);
    printf("Npixel with slope data = %d\n",int(count+0.5));
    printf("<slopeDeviation> (mrad) = %f +- %f\n",sdm,sigma_sdm);
    printf("slopeDeviationMAX (mrad) = %f \n",sdMAX);
    sdm=sumX/count;
    sigma_sdm=(sumX2-sumX*sumX/count)/double(count-1);
    sigma_sdm=sqrt(sigma_sdm);
    printf("<slopeDeviation_X> (mrad) = %f +- %f\n", sdm,sigma_sdm);
    sdm=sumY/count;
    sigma_sdm=(sumY2-sumY*sumY/count)/double(count-1);
    sigma_sdm=sqrt(sigma_sdm);
    printf("<slopeDeviation_Y> (mrad) = %f +- %f\n", sdm,sigma_sdm);
    printf(">>>>> job completed! <<<<<<<\n");
    fflush(stdout);
    file.close();
}


void unitV12(double x1, double y1, double z1, double x2, double y2, double z2){
 // calculates the unit vector P1->P2
 double sum=0.;
 Vservice[0]=x2-x1;
 Vservice[1]=y2-y1;
 Vservice[2]=z2-z1;
 for(int i=0;i<3;i++)
 sum=sum+Vservice[i]*Vservice[i];
 sum=sqrt(sum);
 for(int i=0;i<3;i++)
  Vservice[i]=Vservice[i]/sum;
}


void unitVnIde(double x, double y,double focal){
 //calculates the unit vector normal to the ideal paraboloid surface in x,y
    double z,alpha,r,Dr,vn[3],mvn;
    alpha=atan2(y,x);
    r=sqrt(x*x+y*y);
    z=0.25/focal*r*r;
    Dr=0.5*r/focal;
    vn[0]=-Dr*cos(alpha);
    vn[1]=-Dr*sin(alpha);
    vn[2]=1.0;
    mvn=0.;
    for(int i=0;i<3;i++)
        mvn=mvn+vn[i]*vn[i];
    mvn=sqrt(mvn);
    for(int i=0;i<3;i++)
       Vservice[i]=vn[i]/mvn;
    caronte=z;
}

void pxColor(double val){
    //compute the pixel color to build a color map
    int Blue=0;
    int Green=0;
    int Red=0;
    if(val==0.){
        Blue=255;
        Green=0;
        Red=0;
    }
    else if(val>0. && val< 0.25){
        Blue=255;
        Green=int(255.*4.*val+0.5);
        Red=0;
    }
    else if(val>= 0.25 && val< 0.50){
        Blue=int(255.*(1.-4.*(val-0.25)));
        Green=255;
        Red=0;
    }
    else if(val>= 0.50 && val< 0.75){
        Blue=0;
        Green=255;
        Red=int(255.*4.*(val-0.50));
    }
    else if(val>= 0.75 && val<=1.0){
        Blue=0;
        Green=int(255.*(1.-4.*(val-0.75)));
        Red=255;
    }
    Vservice[0]=Blue;
    Vservice[1]=Green;
    Vservice[2]=Red;
}


inline QImage  cvMatToQImage( const cv::Mat &inMat )
   {
      switch ( inMat.type() )
      {
         // 8-bit, 4 channel
         case CV_8UC4:
         {
            QImage image( inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_RGB32 );

            return image;
         }

         // 8-bit, 3 channel
         case CV_8UC3:
         {
            QImage image( inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_RGB888 );

            return image.rgbSwapped();
         }

         // 8-bit, 1 channel
         case CV_8UC1:
         {
            static QVector<QRgb>  sColorTable;

            // only create our color table once
            if ( sColorTable.isEmpty() )
            {
               for ( int i = 0; i < 256; ++i )
                  sColorTable.push_back( qRgb( i, i, i ) );
            }

            QImage image( inMat.data, inMat.cols, inMat.rows, inMat.step, QImage::Format_Indexed8 );

            image.setColorTable( sColorTable );

            return image;
         }

         default:
            qWarning() << "ASM::cvMatToQImage() - cv::Mat image type not handled in switch:" << inMat.type();
            break;
      }

      return QImage();
   }

   inline QPixmap cvMatToQPixmap( const cv::Mat &inMat )
   {
      return QPixmap::fromImage( cvMatToQImage( inMat ) );
   }
