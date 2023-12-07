/*Author: Marco Montecchi
             ENEA (Italy)
             email: marco.montecchi@enea.it

kFluxMapper controls a GigE camera (by means of the Aravis library) for acquiring
the images needed to the experimental evaluation of the profile of the solar radiation
concentrated by the solar dish in its focal plane where a diffusive plane target is set


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
#include <arv.h>
#include <stdlib.h>
#include <signal.h>
#include <stdio.h>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <fstream>
#include <iostream>
#include <QFileDialog>
#include <QFileInfo>

using namespace std;
using namespace cv;

#include "kfluxmapper.h"
#include "ui_kfluxmapper.h"

//known everywhere
int iSave[3]={0,0,0};
QString PathSave="/run/media/marco/d2quadra/FluxMapper/";
string LogFile="/run/media/marco/d2quadra/FluxMapper/LogFile.txt";
double expTime;
gint payload;
ArvCamera *camera;
ArvStream *stream;
int Width,Height;

kFluxMapper::kFluxMapper(QWidget *parent) : QMainWindow(parent), ui(new Ui::kFluxMapper)
{
    ui->setupUi(this);
    // Connect GUI events with event handlers
    connect(ui->pB_aquire, SIGNAL( clicked() ), this, SLOT( Aquire() ) );
    connect(ui->spinBox_back,SIGNAL(valueChanged()),this,SLOT(Adj_Black()));
    connect(ui->spinBox_oneSun,SIGNAL(valueChanged()),this,SLOT(Adj_OneS()));
    connect(ui->spinBox_flux,SIGNAL(valueChanged()),this,SLOT(Adj_Flux()));
    connect(ui->pB_multiply,SIGNAL(clicked()),this,SLOT(Adj_MulExpTime()));

    //check existing saved imgN.ppm file in PathSave
    QString fext;
    for(int i=0;i<3;i++){
        int iRun=1;
        int iS=1;
        if(i==0)
            fext="Back";
        else if(i==1)
            fext="OneS";
        else if(i==2)
            fext="Flux";
        do{
            QString fileName=PathSave+"img_"+fext+"_"+QString::number(iS)+".ppm";
            QFile fi(fileName);
            if(fi.exists())
                iS++;
            else
                iRun=0;
        }while(iRun==1);
        iS--;
        iSave[i]=iS;
        if(i==0)
            ui->spinBox_back -> setValue(iS);
        else if(i==1)
            ui->spinBox_oneSun -> setValue(iS);
        else if(i==2)
            ui->spinBox_flux -> setValue(iS);
    }

    Width= 2448;
    Height=2050;
    camera = arv_camera_new (NULL); /* Instantiation of the first available camera */
    const char *Cvendor=arv_camera_get_vendor_name (camera);
    const char *Cmodel=arv_camera_get_model_name (camera);
    QString VendorModel=QString(Cvendor)+"\t"+QString(Cmodel);
    ui->lineEdit_camera->setText(VendorModel);
    //set continuous acquisition
    arv_camera_set_acquisition_mode (camera,arv_acquisition_mode_from_string ("ARV_ACQUISITION_MODE_CONTINUOUS"));
    //set software trigger
    arv_camera_set_trigger (camera,"Software");
    //Set mono 14 bit
    arv_camera_set_pixel_format(camera,0x01100025);
    const char *Fpixel=arv_camera_get_pixel_format_as_string (camera);
    ui->lineEdit_depth->setText(Fpixel);
    //set ROI
    arv_camera_set_region (camera, 0, 0, Width, Height);
    /* Set frame rate 5 Hz */
    arv_camera_set_frame_rate (camera, 5.0);
    /* retrieve image payload (number of bytes per image) */
    payload = arv_camera_get_payload (camera);
    ui->lineEdit_width -> setText(QString::number(Width));
    ui->lineEdit_height -> setText(QString::number(Height));

    namedWindow( "Gray image", WINDOW_NORMAL );
    Mat matDisplay(Height,Width,CV_16U);//to avoid warnings
    imshow( "Gray image", matDisplay );

}

kFluxMapper::~kFluxMapper()
{
    delete ui;
}

void kFluxMapper::Adj_Black(){
    iSave[0]=ui->spinBox_back -> value();
}

void kFluxMapper::Adj_OneS(){
    iSave[1]=ui->spinBox_oneSun -> value();
}

void kFluxMapper::Adj_Flux(){
    iSave[2]=ui->spinBox_flux -> value();
}

void kFluxMapper::Aquire(){
    expTime = ui->dSBexpTime -> value();
    //Set Exp. time
    arv_camera_set_exposure_time(camera,expTime);
    /* Create a new stream object */
    stream = arv_camera_create_stream (camera, NULL, NULL);
    if (stream != NULL) {
        /* Push 50 buffer in the stream input buffer queue */
        //for (int i = 0; i < 50; i++)
        arv_stream_push_buffer (stream, arv_buffer_new (payload, NULL));
        /* Start the video stream */
        arv_camera_start_acquisition (camera);
        //send software trigger
        arv_camera_software_trigger (camera);
        ArvBuffer *buffer;
        buffer = arv_stream_pop_buffer(stream);
        if (buffer != NULL) {
            if (arv_buffer_get_status (buffer) == ARV_BUFFER_STATUS_SUCCESS)
                arv_stream_push_buffer (stream, buffer);
            //const unsigned char *imagebuffer = NULL;
            const void *imagebuffer;
            size_t imageSize;
            imagebuffer=arv_buffer_get_data(buffer, &imageSize);
            //salvataggio su matrice openCV
            Mat mat(Height, Width, CV_16U); // mono 14 bit
            Mat matDisplay(Height,Width,CV_16UC3);
            memcpy(mat.data, imagebuffer, Width*Height*(imageSize/Width/Height));
            imwrite(PathSave.toStdString()+"imagine.ppm",mat);
            uint16_t value;
            int valMax=0;
            int valMin=int(pow(2.,14.));
            double magn=pow(2.,16)/pow(2.,14);
            for(int i=0; i<Height; i++){
                for(int j=0; j< Width; j++){
                    value=mat.at<uint16_t>(i,j);
                    valMin=MIN(valMin,value);
                    valMax=MAX(valMax,value);
                    if(value < pow(2.,14)-1){
                        matDisplay.at<Vec3s>(i,j)[0]=int(double(value)*magn);
                        matDisplay.at<Vec3s>(i,j)[1]=int(double(value)*magn);
                        matDisplay.at<Vec3s>(i,j)[2]=int(double(value)*magn);
                    }
                    else{
                        matDisplay.at<Vec3s>(i,j)[0]=0;
                        matDisplay.at<Vec3s>(i,j)[1]=0;
                        matDisplay.at<Vec3s>(i,j)[2]=short(65535);
                    }
                }
            }
            imshow("Gray image",matDisplay);
            waitKey(10);
            ui->lineEdit_MIN -> setText(QString::number(valMin));
            if(valMax < int(pow(2.,14.))-1)
                ui->lineEdit_MAX -> setText(QString::number(valMax));
            else
                ui->lineEdit_MAX -> setText(QString::number(valMax)+"!!!");
            int iType=ui->comboBox -> currentIndex();
            QString fext;
            if(iType==0)
                fext="Back";
            else if(iType==1)
                fext="OneS";
            else if(iType==2)
                fext="Flux";
            Qt::CheckState state1;
            state1 = ui->checkBox -> checkState();
            if( state1 == Qt::Checked ) {
                iSave[iType]++;
                QString fileName = PathSave+"img_"+fext+"_"+QString::number(iSave[iType])+".ppm";
                imwrite(fileName.toStdString(),mat);
                if(iType==0)
                    ui->spinBox_back -> setValue(iSave[iType]);
                else if(iType==1)
                    ui->spinBox_oneSun -> setValue(iSave[iType]);
                else if(iType==2)
                    ui->spinBox_flux -> setValue(iSave[iType]);
                double z= ui->spinBox_z -> value();
                ofstream outfile;
                outfile.open(LogFile, std::ios_base::app);
                outfile << "img_"<<fext.toStdString()<<"_" << iSave[iType] << "\t" << expTime <<"\t"<< z <<endl;
            }
        }
        /* Stop the video stream */
        arv_camera_stop_acquisition (camera);
    }

}


void kFluxMapper::Adj_MulExpTime(){
    expTime = ui->dSBexpTime -> value();
    int factor = ui->dSBexpTimeMulty ->value();
    expTime=expTime*double(factor);
    ui->dSBexpTime -> setValue(expTime);

}
