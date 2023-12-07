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
#include "kfluxmapper.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    kFluxMapper w;
    w.show();

    return a.exec();
}
