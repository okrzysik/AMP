#include "mainwindow.h"

#include <QAction>
#include <QApplication>
#include <QObject>
#include <QtCore>
#include <thread>


int main( int argc, char *argv[] )
{
    // Load the application
    // Q_INIT_RESOURCE(application);
    QApplication app( argc, argv );
    app.setOrganizationName( "" );
    app.setApplicationName( "get_lines" );
    MainWindow mainWin;
#if defined( Q_OS_SYMBIAN )
    mainWin.showMaximized();
#else
    mainWin.show();
#endif
    // Wait for application to finish
    return app.exec();
}
