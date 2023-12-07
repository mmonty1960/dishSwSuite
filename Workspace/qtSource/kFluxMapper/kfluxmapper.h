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
#ifndef KFLUXMAPPER_H
#define KFLUXMAPPER_H

#include <QMainWindow>

namespace Ui {
class kFluxMapper;
}

class kFluxMapper : public QMainWindow
{
    Q_OBJECT

public:
    explicit kFluxMapper(QWidget *parent = 0);
    ~kFluxMapper();

private:
    Ui::kFluxMapper *ui;

public slots:
   void Adj_Black();
   void Adj_OneS();
   void Adj_Flux();
   void Aquire();
   void Adj_MulExpTime();
};

#endif // KFLUXMAPPER_H
