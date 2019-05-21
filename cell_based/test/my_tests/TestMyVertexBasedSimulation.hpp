/*

Copyright (c) 2005-2018, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTMYVERTEXBASEDSIMULATION_HPP_
#define TESTMYVERTEXBASEDSIMULATION_HPP_

/*
 * = Examples showing how to create, run and visualize vertex-based simulations =
 *
 * EMPTYLINE
 *
 * == Introduction ==
 *
 * EMPTYLINE
 *
 * In this tutorial we show how Chaste can be used to create, run and visualize vertex-based simulations.
 * Full details of the mechanical model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 *
 * EMPTYLINE
 *
 * == The test ==
 *
 * EMPTYLINE
 *
 * As in previous cell-based Chaste tutorials, we begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

/* The remaining header files define classes that will be used in the cell-based
 * simulation. We have encountered some of these header files in previous cell-based
 * Chaste tutorials. */
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the cell cycle model. */
#include "UniformG1GenerationalCellCycleModel.hpp"
/* The next two header files define a helper class for generating suitable meshes: one planar and one periodic. */
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
/* The next header file defines a vertex-based {{{CellPopulation}}} class.*/
#include "VertexBasedCellPopulation.hpp"
/* The next header file defines a force law for describing the mechanical interactions
 * between neighbouring cells in the cell population, subject to each vertex.
 */
#include "NagaiHondaForce.hpp"
/* This force law assumes that cells possess a "target area" property which determines the size of each
 * cell in the simulation. In order to assign target areas to cells and update them in each time step, we need
 * the next header file.
 */
#include "SimpleTargetAreaModifier.hpp"
/* The next header file defines a boundary condition for the cells.*/
#include "PlaneBoundaryCondition.hpp"
/* The next header file defines a cell killer, which specifies how cells are removed from the simulation.*/
#include "PlaneBasedCellKiller.hpp"

/* Finally, we include a header that enforces running this test only on one process. */
#include "FakePetscSetup.hpp"


class TestRunningVertexBasedSimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    
    void TestMonolayer()
    {
        
        HoneycombVertexMeshGenerator generator(3, 3);    // Parameters are: cells across, cells up
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
      
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MyVertexBasedSimulation");
        simulator.SetEndTime(1.0);
        
        simulator.SetSamplingTimestepMultiple(50);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);
       
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();
    }

};

#endif /* TESTMYVERTEXBASEDSIMULATION_HPP_ */
