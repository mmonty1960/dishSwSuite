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
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Widget w;
    w.show();

    return a.exec();
}
