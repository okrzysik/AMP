#include <QAction>
#include <QCloseEvent>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QImage>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QSqlDatabase>
#include <QSqlTableModel>
#include <QStatusBar>
#include <QtGui>

#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <thread>
#include <vector>

#include "mainwindow.h"


extern "C" {
#include <cassert>
}


// Function to add an action (with connection) to a menu
#define ADD_MENU_ACTION( menu, string, arg )                                   \
    do {                                                                       \
        QAction *action = new QAction( string, this );                         \
        connect( action, SIGNAL( triggered() ), signalMapper, SLOT( map() ) ); \
        signalMapper->setMapping( action, arg );                               \
        menu->addAction( action );                                             \
    } while ( 0 )


#define ASSERT( EXP )                                                                     \
    do {                                                                                  \
        if ( !( EXP ) ) {                                                                 \
            std::stringstream stream;                                                     \
            stream << "Failed assertion: " << #EXP << " " << __FILE__ << " " << __LINE__; \
            throw std::logic_error( stream.str() );                                       \
        }                                                                                 \
    } while ( 0 )


// Class to draw a box
class DrawBoxClass : public QWidget
{
public:
    DrawBoxClass( QWidget *parent, int boarder = 0 ) : QWidget( parent ), d_boarder( boarder ) {}
    // Set first coordinate
    void set_p1( QPoint pos )
    {
        auto p = mapFromGlobal( pos );
        auto s = parentWidget()->size();
        int x  = std::min( std::max( p.x(), 0 ), s.width() );
        int y  = std::min( std::max( p.y(), 0 ), s.height() );
        d_p1   = QPoint( x, y );
    }
    void set_p2( QPoint pos )
    {
        auto p = mapFromGlobal( pos );
        auto s = parentWidget()->size();
        int x  = std::min( std::max( p.x(), 0 ), s.width() );
        int y  = std::min( std::max( p.y(), 0 ), s.height() );
        d_p2   = QPoint( x, y );
    }
    // Update
    void redraw()
    {
        const QRect &geom = parentWidget()->geometry();
        setGeometry( QRect( 0, 0, geom.width(), geom.height() ) );
        this->update();
    }
    // Draw the rectangle
    void paintEvent( QPaintEvent * ) override
    {
        QPainter p( this );
        QPen pen( QColor( 60, 0, 0 ) );
        pen.setWidth( 3 );
        p.setPen( pen );
        p.setRenderHint( QPainter::Antialiasing );
        int x = std::min( d_p1.x(), d_p2.x() );
        int y = std::min( d_p1.y(), d_p2.y() );
        int w = std::max( d_p1.x(), d_p2.x() ) - x;
        int h = std::max( d_p1.y(), d_p2.y() ) - y;
        p.drawRect( QRect( x, y, w, h ) );
    }

protected:
    int d_boarder;
    QPoint d_p1, d_p2;
};


/***********************************************************************
 * Constructor/destructor                                               *
 ***********************************************************************/
MainWindow::MainWindow()
    : ThreadedSlotsClass(),
      mainMenu( nullptr ),
      backButton( nullptr ),
      processorButton( nullptr ),
      imageWindow( nullptr )
{
    QWidget::setWindowTitle( QString( "GetLines" ) );

    // Back button
    backButton = new QPushButton( "Back" );
    connect( backButton, SIGNAL( released() ), this, SLOT( backButtonPressed() ) );
    backButton->hide();

    // Create the image plot
    imageWindow = new QLabel;
    // imageWindow->setMinimumSize( 100, 100 );
    // imageWindow->setSizePolicy(QSizePolicy::Maximum, QSizePolicy::Maximum);

    auto *layout = new QVBoxLayout;
    layout->setMargin( 0 );
    layout->setContentsMargins( QMargins( 0, 0, 0, 0 ) );
    layout->setSpacing( 0 );
    layout->addWidget( imageWindow );
    setCentralWidget( new QWidget );
    centralWidget()->setLayout( layout );
    gbox.reset( new DrawBoxClass( imageWindow, 2 ) );
    gbox->setVisible( false );

    createActions();
    createMenus();
    createToolBars();
    createStatusBar();

    readSettings();

    // Create resize event
    resizeTimer.setSingleShot( true );
    connect( &resizeTimer, SIGNAL( timeout() ), SLOT( resizeDone() ) );

    setUnifiedTitleAndToolBarOnMac( true );
    reset();
    updateDisplay();
    qApp->processEvents();
}
MainWindow::~MainWindow()
{
    // Close the file
    close();
    // Disconnect signals created by createActions
    disconnect( openAct, SIGNAL( triggered() ), this, SLOT( open() ) );
    disconnect( closeAct, SIGNAL( triggered() ), this, SLOT( close() ) );
    disconnect( resetAct, SIGNAL( triggered() ), this, SLOT( reset() ) );
    disconnect( exclusiveAct, SIGNAL( triggered() ), this, SLOT( exclusiveFun() ) );
    disconnect( subfunctionsAct, SIGNAL( triggered() ), this, SLOT( subfunctionFun() ) );
    disconnect( exitAct, SIGNAL( triggered() ), this, SLOT( exit() ) );
    disconnect( aboutAct, SIGNAL( triggered() ), this, SLOT( about() ) );
    disconnect( savePerformanceTimers, SIGNAL( triggered() ), this, SLOT( savePerformance() ) );
    disconnect( runUnitTestAction, SIGNAL( triggered() ), this, SLOT( runUnitTestsSlot() ) );
    // Disconnect signals created by createToolbars
    disconnect( exclusiveButton, SIGNAL( released() ), this, SLOT( exclusiveFun() ) );
    disconnect( subfunctionButton, SIGNAL( released() ), this, SLOT( subfunctionFun() ) );
    // Delete objects
    delete mainMenu;
}
void MainWindow::closeEvent( QCloseEvent *event )
{
    writeSettings();
    close();
    event->accept();
}
void MainWindow::close()
{
    /*d_data = TimerMemoryResults();
    d_dataTimer.clear();
    d_dataTrace.clear();
    traceWindow.reset();
    memoryWindow.reset();*/
    QWidget::setWindowTitle( QString( "GetLines" ) );
    reset();
}
void MainWindow::exit() { qApp->quit(); }


/***********************************************************************
 * Reset the view                                                       *
 ***********************************************************************/
void MainWindow::reset() { updateDisplay(); }


/***********************************************************************
 * Help                                                                 *
 ***********************************************************************/
void MainWindow::about()
{
    QMessageBox::about(
        this, tr( "About Application" ), tr( "This is a program to load lines from an image" ) );
}


/***********************************************************************
 * Save the timers for the application                                  *
 ***********************************************************************/
void MainWindow::savePerformance()
{
    /*std::string filename = QFileDialog::getSaveFileName(
        this, "Name of file to save", lastPath.c_str(), "*.0.timer *.1.timer" )
                               .toStdString();
    if ( !filename.empty() ) {
        size_t i = filename.rfind( '.', filename.rfind( ".timer" ) - 1 );
        filename = filename.substr( 0, i );
        PROFILE_SAVE( filename );
    }*/
}


/***********************************************************************
 * Load timer data                                                      *
 ***********************************************************************/
void MainWindow::open()
{
    QString filename = QFileDialog::getOpenFileName( this, "Select the image", lastPath.c_str() );
    loadFile( filename.toStdString() );
}
void MainWindow::loadFile( const std::string &filename, bool showFailure )
{
    if ( !filename.empty() ) {
        // Load the image
        QImage image( QString( filename.data() ) );
        data.resize( image.width(), image.height() );
        data.fill( 0 );
        for ( int i = 0; i < image.width(); i++ ) {
            for ( int j = 0; j < image.height(); j++ ) {
                auto rgb     = image.pixel( i, j );
                data( i, j ) = AMP::RGBA( qRed( rgb ), qGreen( rgb ), qBlue( rgb ) );
            }
        }
        // Set the title window
        QWidget::setWindowTitle( filename.data() );
        // Update the display
        updateDisplay();
        auto resizeEvent = new QResizeEvent( size(), size() );
        QCoreApplication::postEvent( this, resizeEvent );
    }
}


/***********************************************************************
 * Update the display                                                   *
 ***********************************************************************/
template<class TYPE>
inline TYPE getTableData( const std::vector<TYPE> &x, int rank )
{
    TYPE y;
    if ( rank == -1 ) {
        y = static_cast<TYPE>( mean( x ) );
    } else if ( rank == -2 ) {
        y = static_cast<TYPE>( min( x ) );
    } else if ( rank == -3 ) {
        y = static_cast<TYPE>( max( x ) );
    } else {
        y = x[rank];
    }
    return y;
}
void MainWindow::updateDisplay()
{
    if ( data.empty() ) {
        // No data, hide objects
        imageWindow->setVisible( false );
        return;
    }
    // Plot the image
    QImage image( (uchar *) data.data(), data.size( 0 ), data.size( 1 ), QImage::Format_RGB32 );
    QPixmap imagePixelMap( QPixmap::fromImage( image ) );
    int w       = imageWindow->rect().width();
    int h       = imageWindow->rect().height();
    auto image2 = imagePixelMap.scaled( w, h, Qt::KeepAspectRatio, Qt::SmoothTransformation );
    imageWindow->setPixmap( image2 );
    imageWindow->setVisible( true );
}


/***********************************************************************
 * Back button                                                          *
 ***********************************************************************/
void MainWindow::backButtonPressed()
{
    /*if ( callList.empty() )
        return;
    callList.pop_back();
    updateDisplay();*/
}


/***********************************************************************
 * Create the actions                                                   *
 ***********************************************************************/
void MainWindow::createActions()
{
    openAct = new QAction( QIcon( ":/images/open.png" ), tr( "&Open..." ), this );
    openAct->setShortcuts( QKeySequence::Open );
    openAct->setStatusTip( tr( "Open an existing file" ) );
    connect( openAct, SIGNAL( triggered() ), this, SLOT( open() ) );

    closeAct = new QAction( QIcon( ":/images/new.png" ), tr( "&Close" ), this );
    closeAct->setShortcuts( QKeySequence::Close );
    closeAct->setStatusTip( tr( "Create a new file" ) );
    connect( closeAct, SIGNAL( triggered() ), this, SLOT( close() ) );

    resetAct = new QAction( QIcon( ":/images/refresh.png" ), tr( "&Reset" ), this );
    resetAct->setShortcuts( QKeySequence::Refresh );
    resetAct->setStatusTip( tr( "Reset the view" ) );
    connect( resetAct, SIGNAL( triggered() ), this, SLOT( reset() ) );

    exitAct = new QAction( tr( "E&xit" ), this );
    exitAct->setShortcuts( QKeySequence::Quit );
    exitAct->setStatusTip( tr( "Exit the application" ) );
    connect( exitAct, SIGNAL( triggered() ), this, SLOT( exit() ) );

    aboutAct = new QAction( tr( "&About" ), this );
    aboutAct->setStatusTip( tr( "Show the application's About box" ) );
    connect( aboutAct, SIGNAL( triggered() ), this, SLOT( about() ) );

    savePerformanceTimers = new QAction( tr( "&Save performance" ), this );
    savePerformanceTimers->setStatusTip( tr( "Save performance timers" ) );
    connect( savePerformanceTimers, SIGNAL( triggered() ), this, SLOT( savePerformance() ) );

    runUnitTestAction = new QAction( tr( "&Run unit tests" ), this );
    runUnitTestAction->setStatusTip( tr( "Run unit tests" ) );
    connect( runUnitTestAction, SIGNAL( triggered() ), this, SLOT( runUnitTestsSlot() ) );
}

void MainWindow::createMenus()
{
    mainMenu = new QMenuBar( nullptr );
    fileMenu = mainMenu->addMenu( tr( "&File" ) );
    fileMenu->addAction( openAct );
    fileMenu->addAction( closeAct );
    fileMenu->addAction( resetAct );
    // fileMenu->addAction(saveAct);
    // fileMenu->addAction(saveAsAct);
    // fileMenu->addSeparator();
    fileMenu->addAction( exitAct );

    // editMenu = mainMenu->addMenu(tr("&Edit"));
    // editMenu->addAction(cutAct);
    // editMenu->addAction(copyAct);
    // editMenu->addAction(pasteAct);

    mainMenu->addSeparator();

    helpMenu = mainMenu->addMenu( tr( "&Help" ) );
    helpMenu->addAction( aboutAct );

    // Create developer menu on rhs
    auto *developer_bar = new QMenuBar( this );
    QMenu *developer    = developer_bar->addMenu( tr( "&Developer" ) );
    developer->addAction( savePerformanceTimers );
    developer->addAction( runUnitTestAction );
    mainMenu->setCornerWidget( developer_bar );

// In Ubuntu 14.04 with qt5 the window's menu bar goes missing entirely
// if the user is running any desktop environment other than Unity
// (in which the faux single-menubar appears). The user has a
// workaround, to remove the appmenu-qt5 package, but that is
// awkward and the problem is so severe that it merits disabling
// the system menubar integration altogether. Like this:
#if defined( Q_OS_LINUX ) && QT_VERSION >= 0x050000
    mainMenu->setNativeMenuBar( false );      // fix #1039
    developer_bar->setNativeMenuBar( false ); // fix #1039
#endif

    setMenuBar( mainMenu );
}

void MainWindow::createToolBars()
{
    fileToolBar = addToolBar( tr( "File" ) );
    fileToolBar->addAction( closeAct );
    fileToolBar->addAction( openAct );
    fileToolBar->addSeparator();
    fileToolBar->addAction( resetAct );
}

void MainWindow::createStatusBar() { statusBar()->showMessage( tr( "Ready" ) ); }

void MainWindow::readSettings()
{
    QSettings settings( "AMP", "get_image" );
    QPoint pos = settings.value( "pos", QPoint( 200, 200 ) ).toPoint();
    QSize size = settings.value( "size", QSize( 400, 400 ) ).toSize();
    resize( size );
    move( pos );
}

void MainWindow::writeSettings()
{
    QSettings settings( "AMP", "get_image" );
    settings.setValue( "pos", pos() );
    settings.setValue( "size", size() );
}


/***********************************************************************
 * Resize functions                                                     *
 ***********************************************************************/
void MainWindow::resizeEvent( QResizeEvent *e )
{
    resizeTimer.start( 200 );
    QMainWindow::resizeEvent( e );
}
void MainWindow::resizeDone() { updateDisplay(); }


/***********************************************************************
 * Mouse events                                                         *
 ***********************************************************************/
inline std::ostream &operator<<( std::ostream &out, const QPoint &p )
{
    out << "(" << p.x() << "," << p.y() << ")";
    return out;
}
void MainWindow::mousePressEvent( QMouseEvent *event )
{
    gbox->set_p1( mapToGlobal( event->pos() ) );
    gbox->set_p2( mapToGlobal( event->pos() ) );
    gbox->redraw();
    gbox->setVisible( true );
}
void MainWindow::mouseMoveEvent( QMouseEvent *event )
{
    gbox->set_p2( mapToGlobal( event->pos() ) );
    gbox->redraw();
}
void MainWindow::mouseReleaseEvent( QMouseEvent *event )
{
    gbox->set_p2( mapToGlobal( event->pos() ) );
    gbox->setVisible( false );
}


/***********************************************************************
 * Run the unit tests                                                   *
 ***********************************************************************/
void MainWindow::runUnitTestsSlot()
{
    // Get the filename to test
    QString filename = QFileDialog::getOpenFileName(
        this, "Select the timer file to test", lastPath.c_str(), "*.0.timer *.1.timer" );
    // Run the unit tests
    bool pass = runUnitTests( filename.toStdString() );
    if ( pass ) {
        QMessageBox::information( this, tr( "All tests passed" ), tr( "All tests passed" ) );
    } else {
        QMessageBox::information( this, tr( "Some tests failed" ), tr( "Some tests failed" ) );
    }
}
void MainWindow::callLoadFile()
{
    loadFile( unitTestFilename, false );
    qApp->processEvents();
}
bool MainWindow::runUnitTests( const std::string &filename )
{
    /*    // Try to load the file
        unitTestFilename = filename;
        callFun( std::bind( &MainWindow::callLoadFile, this ) );
        if ( d_dataTimer.empty() ) {
            printf( "   Failed to load file %s\n", filename.c_str() );
            return false;
        }
        // Select the first cell
        callFun( std::bind( &MainWindow::callSelectCell, this ) );

        // Run the trace unit tests
        if ( hasTraceData() ) {
            callFun( std::bind( &MainWindow::traceFun, this ) );
            if ( !traceWindow->runUnitTests() )
                return false;
            callFun( std::bind( &MainWindow::closeTrace, this ) );
        }
        // Try to close the file
        callFun( std::bind( &MainWindow::close, this ) );
        if ( !d_dataTimer.empty() ) {
            printf( "   Failed call to close\n" );
            return false;
        }
        return true;*/
    return true;
}
