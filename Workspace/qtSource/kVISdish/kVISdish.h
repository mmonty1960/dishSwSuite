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
#ifndef KVISDISH_H
#define KVISDISH_H

#include <QWidget>
#include <iostream>
#include <fstream>
#include <QTextStream>
#include <QFile>

/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia
Roma - Italy*/


namespace Ui {
class Widget;
}

class Widget : public QWidget
{
    Q_OBJECT

public:
    explicit Widget(QWidget *parent = 0);
    ~Widget();
    void kPE();

public slots:
    void plotBlob();
    void takeOnePicture();
    void Process();
private slots:
    void on_pushButton_clicked();
    void on_pushButton_slaveFolder_clicked();


private:
    Ui::Widget *ui;

};

#endif // KVISDISH_H
