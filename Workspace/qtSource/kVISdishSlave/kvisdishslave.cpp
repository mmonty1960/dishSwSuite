/*Author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it

kVISdishisSlave is the slave software for infield experimental evaluation of 3D-shape of
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
#include "kvisdishslave.h"
#include "ui_kvisdishslave.h"
#include <QApplication>
#include <QtGui>
#include <QStringList>
#include <QFileSystemWatcher>
#include <unistd.h>
#include <QMessageBox>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

using namespace std;
using namespace cv;

QString Path = "/home/marco/Immagini/blobPattern";
QString fileCTRL = Path+ "/currentBlob.ctrl";

int iOkFS=0;


kVISdishSlave::kVISdishSlave(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::kVISdishSlave)
{
    ui->setupUi(this);

    qfsw = new QFileSystemWatcher(this);
    qfsw -> addPath(fileCTRL);

    ui->lineEdit_image->setText(Path);

    connect(qfsw, SIGNAL(fileChanged(QString)),this, SLOT( fC(QString)));

    namedWindow("oneBlob",WINDOW_NORMAL);

    fC(Path);
}

kVISdishSlave::~kVISdishSlave()
{
    delete ui;
}


void kVISdishSlave::fC(const QString & ){
    int iB1L2,iBW,rBlob,hight,width;
    Point2f Pt1,Pt2;
    cv::Scalar scalar;
    Mat mat;

    QFile filectrl(fileCTRL);
    while (!filectrl.open(QIODevice::ReadOnly | QIODevice::Text))
      usleep(10);
    QTextStream stream ( &filectrl );
    QString info = stream.readLine();
    printf("info: %s\n",(info.toStdString()).c_str());
    if(info.contains("fullscreen", Qt::CaseInsensitive) || info.contains("normalscreen", Qt::CaseInsensitive)){
        if(info.contains("fullscreen", Qt::CaseInsensitive)){
            if(iOkFS==0) setWindowProperty("oneBlob", WND_PROP_FULLSCREEN, WINDOW_FULLSCREEN);
            iOkFS=1;
        }
        if(info.contains("normalscreen", Qt::CaseInsensitive)){
            if(iOkFS==1) {
                destroyWindow("oneBlob");
                namedWindow("oneBlob",WINDOW_NORMAL);
            }
            iOkFS=0;
        }
        stream >> hight >> width;
        printf("hight = %d  width = %d\n",hight,width);
        while (!stream.atEnd() && hight >0 && width >0){
            stream >> iB1L2 >> iBW >> rBlob;
            printf("iB1L2 = %d  iBW = %d  rBlob = %d\n",iB1L2,iBW,rBlob);
            if(iB1L2 >0){
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
                if(iB1L2==1 && rBlob >0){
                    stream >> Pt1.x >> Pt1.y;
                    printf("Blob: Pt1.x = %f  Pt1.y =%f\n",Pt1.x,Pt1.y);
                    circle(mat, Pt1, rBlob, scalar, -1);
                } else if(iB1L2==2 && rBlob >0){
                    stream >> Pt1.x >> Pt1.y >>Pt2.x >> Pt2.y;
                    printf("Line: Pt1.x = %f  Pt1.y =%f Pt2.x = %f  Pt2.y =%f\n",Pt1.x,Pt1.y,Pt2.x,Pt2.y);
                    line(mat,Pt1,Pt2,scalar,2*rBlob,8,0);
                }
                imshow("oneBlob",mat);
                waitKey(10);
            }
        }
    }
    fflush(stdout);
    filectrl.close();
    mat.release();
//    QFileInfo checkFile(fileCTRL);
//    while(!checkFile.exists())
//        usleep(10);
//    qfsw->addPath(fileCTRL);
}
