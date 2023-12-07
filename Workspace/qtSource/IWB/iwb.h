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
#ifndef IWB_H
#define IWB_H

//#include <QWidget>
//#include "opencv2/core/core.hpp"
//#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/imgproc/imgproc.hpp"
//#include "opencv2/calib3d/calib3d.hpp"
#include "ui_iwb.h"

namespace Ui {
class ImageWorkBench;
}

class ImageWorkBench : public QWidget
{
    Q_OBJECT

public:
    explicit ImageWorkBench(QWidget *parent = 0);
    ~ImageWorkBench();

private:
    Ui::ImageWorkBench *ui;

public slots:
    void selectDir();
    void compute();
    void selectImg();
    void viewImg(QString);
    void setIntensify();
    void readGL();
    void saveBG();
    void saveOS();
    void saveFX();
    void setROI();
    void build();
};

#endif // IWB_H
