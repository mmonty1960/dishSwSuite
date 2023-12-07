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
#include <QApplication>
#include "SimulDish.h"
 
int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    SimulDish *dialog = new SimulDish;
    dialog->show();
    return app.exec();
}
