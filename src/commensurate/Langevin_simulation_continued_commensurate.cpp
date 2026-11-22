//
//  main.cpp
//  Langevin dynmaics for XXZ spin chain with dissipation
//
//  Created by smajumdar on 26/03/2021.
//

#include <langevin/commensurate.hpp>
#include <langevin/cli.hpp>

int main(int argc, char* argv[])
{
    return langevin::cli::run<true>(argc, argv, true);
}
