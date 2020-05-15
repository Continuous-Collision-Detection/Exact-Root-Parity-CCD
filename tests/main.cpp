////////////////////////////////////////////////////////////////////////////////
// Keep this file empty, and implement unit tests in separate compilation units!
////////////////////////////////////////////////////////////////////////////////

// Catch2 Documentation: https://github.com/catchorg/Catch2/tree/master/docs

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>


int main(int argc, char* argv[])
{
    Catch::Session session; // There must be exactly one instance

    // Build a new parser on top of Catch's
    using namespace Catch::clara;
    auto cli = session.cli();
    session.cli(cli);

    int returnCode = session.applyCommandLine(argc, argv);
    if (returnCode != 0) // Indicates a command line error
        return returnCode;

    return session.run();
}
