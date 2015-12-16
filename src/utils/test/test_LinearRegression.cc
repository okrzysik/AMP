#include "utils/LinearRegression.h"
#include "utils/AMPManager.h"
#include "utils/UnitTest.h"

void linearRegressionTest( AMP::UnitTest *ut )
{
    double xArray[] = { 71., 73., 64., 65., 61., 70., 65., 72., 63., 67., 64. };
    double yArray[] = { 160., 183., 154., 168., 159., 180., 145., 210., 132., 168., 141. };
    std::vector<double> xVector( xArray, xArray + sizeof( xArray ) / sizeof( double ) );
    std::vector<double> yVector( yArray, yArray + sizeof( yArray ) / sizeof( double ) );

    // std::cout<<"Linear Regression Test"<<std::endl;
    AMP::LinearRegression linReg( xVector, yVector );
    double CoefDeterm  = linReg.getCoefDeterm();
    double CoefCorrel  = linReg.getCoefCorrel();
    double StdErrorEst = linReg.getStdErrorEst();
    double A           = linReg.getA();
    double B           = linReg.getB();
    // std::cout<<"Number of x,y pairs = "<<linReg.items()<<std::endl;
    // std::cout<<linReg<<std::endl;
    // std::cout<<"Coefficient of Determination = "<<linReg.getCoefDeterm()<<std::endl;
    // std::cout<<"Coefficient of Correlation = "<<linReg.getCoefCorrel()<<std::endl;
    // std::cout<<"Standard Error of Estimate = "<<linReg.getStdErrorEst()<<std::endl;

    // Matlab answer
    // >> xx = [  71,  73,  64,  65,  61,  70,  65,  72,  63,  67,  64 ];
    // >> yy = [ 160, 183, 154, 168, 159, 180, 145, 210, 132, 168, 141 ];
    // >> BA = polyfit(xx,yy,1);
    // >> A = BA(2)
    // >> B = BA(1)
    // >> R = corrcoef(xx,yy);
    // >> coefD = R(1,2)^2
    // >> coefC = R(1,2)
    // >> err = yy - polyval(BA, xx);
    // >> stdError = std(err)
    double MatlabCoefDeterm  = 0.556260166947567;
    double MatlabCoefCorrel  = 0.745828510414805;
    double MatlabStdErrorEst = 14.6225187828758;
    double MatlabA           = -106.791666666665;
    double MatlabB           = 4.0472222222222;
    // std::cout<<"A - MatlabA = "<<A - MatlabA<<"\n";
    // std::cout<<"B - MatlabB = "<<B - MatlabB<<"\n";
    // std::cout<<"CoefDeterm - MatlabCoefDeterm = "<<CoefDeterm - MatlabCoefDeterm<<"\n";
    // std::cout<<"CoefDeterm - MatlabCoefCorrel = "<<CoefCorrel - MatlabCoefCorrel<<"\n";
    // std::cout<<"StdErrorEst - MatlabStdErrorEst = "<<StdErrorEst - MatlabStdErrorEst<<"\n";
    std::string thisTest( "testLinearRegression" );
    double CUSTOM_EPSILON = 1.0e-10;
    if ( fabs( A - MatlabA ) > CUSTOM_EPSILON || fabs( B - MatlabB ) > CUSTOM_EPSILON ||
         fabs( CoefDeterm - MatlabCoefDeterm ) > CUSTOM_EPSILON ||
         fabs( CoefCorrel - MatlabCoefCorrel ) > CUSTOM_EPSILON ||
         fabs( StdErrorEst - MatlabStdErrorEst ) > CUSTOM_EPSILON )
        ut->failure( thisTest );
    else
        ut->passes( thisTest );
}

int main( int argc, char *argv[] )
{
    AMP::AMPManager::startup( argc, argv );
    AMP::UnitTest ut;

    try {
        linearRegressionTest( &ut );
    }
    catch ( std::exception &err ) {
        std::cout << "ERROR: While testing " << argv[0] << err.what() << std::endl;
        ut.failure( "ERROR: While testing" );
    }
    catch ( ... ) {
        std::cout << "ERROR: While testing " << argv[0] << "An unknown exception was thrown."
                  << std::endl;
        ut.failure( "ERROR: While testing" );
    }

    ut.report();

    int num_failed = ut.NumFailGlobal();
    AMP::AMPManager::shutdown();
    return num_failed;
}
