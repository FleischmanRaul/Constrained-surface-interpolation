#include "MainWindow.h"
#include <iostream>

using namespace std;

namespace cagd
{
    MainWindow::MainWindow(QWidget *parent): QMainWindow(parent)
    {
        setupUi(this);

    /*

      the structure of the main window's central widget

     *---------------------------------------------------*
     |                 central widget                    |
     |                                                   |
     |  *---------------------------*-----------------*  |
     |  |     rendering context     |   scroll area   |  |
     |  |       OpenGL widget       | *-------------* |  |
     |  |                           | | side widget | |  |
     |  |                           | |             | |  |
     |  |                           | |             | |  |
     |  |                           | *-------------* |  |
     |  *---------------------------*-----------------*  |
     |                                                   |
     *---------------------------------------------------*

    */
        _side_widget = new SideWidget(this);

        _scroll_area = new QScrollArea(this);
        _scroll_area->setWidget(_side_widget);
        _scroll_area->setSizePolicy(_side_widget->sizePolicy());
        _scroll_area->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

        _gl_widget = new GLWidget(this);

        centralWidget()->setLayout(new QHBoxLayout());
        centralWidget()->layout()->addWidget(_gl_widget);
        centralWidget()->layout()->addWidget(_scroll_area);

        // creating a signal slot mechanism between the rendering context and the side widget
        connect(_side_widget->rotate_x_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_angle_x(int)));
        connect(_side_widget->rotate_y_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_angle_y(int)));
        connect(_side_widget->rotate_z_slider, SIGNAL(valueChanged(int)), _gl_widget, SLOT(set_angle_z(int)));

        connect(_side_widget->zoom_factor_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_zoom_factor(double)));

        connect(_side_widget->trans_x_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_trans_x(double)));
        connect(_side_widget->trans_y_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_trans_y(double)));
        connect(_side_widget->trans_z_spin_box, SIGNAL(valueChanged(double)), _gl_widget, SLOT(set_trans_z(double)));

        connect(_side_widget->curveEvolveCheckBox,  SIGNAL(clicked(bool)),       _gl_widget, SLOT(set_evolve_state(bool)));
        connect(_side_widget->surfaceEvolveCheckBox,SIGNAL(clicked(bool)),       _gl_widget, SLOT(set_evolve_state2(bool)));
        connect(_side_widget->curveEnergyType,      SIGNAL(activated(int)),      _gl_widget, SLOT(set_curve_energy_type(int)));
        connect(_side_widget->curveType,            SIGNAL(activated(int)),      _gl_widget, SLOT(set_curve_type(int)));
        connect(_side_widget->createCurvePopulation,SIGNAL(pressed()),           _gl_widget, SLOT(create_curve_population()));
        connect(_side_widget->populationCount,      SIGNAL(valueChanged(int)),   _gl_widget, SLOT(set_population_count(int)));
        connect(_side_widget->maxMaturityLevel,     SIGNAL(valueChanged(int)),   _gl_widget, SLOT(set_max_maturity_level(int)));
        connect(_side_widget->matingPoolCount,      SIGNAL(valueChanged(int)),   _gl_widget, SLOT(set_mating_pool_size(int)));
        connect(_side_widget->threshold,            SIGNAL(valueChanged(double)),_gl_widget, SLOT(set_threshold(double)));

        connect(_side_widget->createCurvePopulation,SIGNAL(pressed()),           this, SLOT(enableCurveEvolveCheckbox()));

        connect(_gl_widget, SIGNAL(valueChanged(double)), _side_widget->epsDoubleSpinBox, SLOT(setValue(double)));
    }

    //--------------------------------
    // implementation of private slots
    //--------------------------------
    void MainWindow::on_action_Quit_triggered()
    {
        qApp->exit(0);
    }

    void MainWindow::enableCurveEvolveCheckbox()
    {
        _side_widget->curveEvolveCheckBox->setEnabled(true);
    }
}
