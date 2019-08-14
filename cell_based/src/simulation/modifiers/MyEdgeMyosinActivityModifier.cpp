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

#include "MyEdgeMyosinActivityModifier.hpp"
#include "VertexElement.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <math.h>

MyEdgeMyosinActivityModifier::MyEdgeMyosinActivityModifier()
    : AbstractCellBasedSimulationModifier<2>(),
      mReferenceTargetArea(1.0)
{
}

MyEdgeMyosinActivityModifier::~MyEdgeMyosinActivityModifier()
{
}

void MyEdgeMyosinActivityModifier::UpdateAtEndOfTimeStep(AbstractCellPopulation<2, 2>& rCellPopulation)
{
    UpdateEdgeMyosinActivitiesAndSurfaceAreaHistories(rCellPopulation);
}

void MyEdgeMyosinActivityModifier::SetupSolve(AbstractCellPopulation<2, 2>& rCellPopulation, std::string outputDirectory)
{
}


double MyEdgeMyosinActivityModifier::GetReferenceTargetArea()
{
    return mReferenceTargetArea;
}

void MyEdgeMyosinActivityModifier::SetReferenceTargetArea(double referenceTargetArea)
{
    assert(referenceTargetArea >= 0.0);
    mReferenceTargetArea = referenceTargetArea;
}

void MyEdgeMyosinActivityModifier::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<ReferenceTargetArea>" << mReferenceTargetArea << "</ReferenceTargetArea>\n";

    // Next, call method on direct parent class
    AbstractCellBasedSimulationModifier<2>::OutputSimulationModifierParameters(rParamsFile);
}

// My changes
void MyEdgeMyosinActivityModifier::UpdateEdgeMyosinActivitiesAndSurfaceAreaHistories(AbstractCellPopulation<2,2>& rCellPopulation)
{
    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*> (&rCellPopulation);
    for (typename VertexMesh<2,2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
            elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
            ++ elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        double previous_element_primeter = elem_iter->GetHistoricSurfaceArea();
        CellPtr p_cell= p_cell_population->GetCellUsingLocationIndex(elem_index);

        // Update MyosinActivity
        double current_time = SimulationTime::Instance()->GetTime();

        double myosin_activity = 0.0;
        double new_myosin_activity = 0.0;
        double edge_myosin_activity = 0.0;
        double new_edge_myosin_activity = 0.0;

        // modify Cell Myosin Activity
        if (current_time<(15.0- 1e-10))
        {
            myosin_activity = p_cell->GetMyosinActivity();
            new_myosin_activity = myosin_activity+ SimulationTime::Instance()->GetTimeStep()*(2.0*pow(previous_element_primeter,8.0)/(1.0*pow(2.063,8.0)+pow(previous_element_primeter,8.0))-myosin_activity);// notice 'pow(2.063,8.0)' in the statement
            if (current_time < 15.0 -1e-10)
                new_myosin_activity = 1.0;
            p_cell->SetMyosinActivity(new_myosin_activity);
        }
        // initialise Edge Myosin Activity
        else if (current_time < (15.0 + 1e-10))
        {
            myosin_activity = p_cell->GetMyosinActivity();
            new_myosin_activity = myosin_activity+ SimulationTime::Instance()->GetTimeStep()*(2.0*pow(previous_element_primeter,8.0)/(1.0*pow(2.063,8.0)+pow(previous_element_primeter,8.0))-myosin_activity);
            elem_iter->InitializeEdgeMyosinActivities(new_myosin_activity);
        }
        // modify Edge Myosin Activity
        else
        {
            for (unsigned i = 0; i < elem_iter->GetNumNodes(); ++i)
            {
                edge_myosin_activity = elem_iter->GetEdgeMyosinActivity(i);
                double previous_edge_length = elem_iter->GetHistoricEdgeLength(i);
                new_edge_myosin_activity = edge_myosin_activity+ SimulationTime::Instance()->GetTimeStep()*(2.0*pow(6*previous_edge_length,8.0)/(pow(2.063,8.0)+pow(6*previous_edge_length,8.0))-edge_myosin_activity);
                //double previous_edge_length = 0.6204;
                //new_edge_myosin_activity = edge_myosin_activity+ SimulationTime::Instance()->GetTimeStep()*(2.0*pow(previous_element_primeter,8.0)/(pow(2.063,8.0)+pow(previous_element_primeter,8.0))-edge_myosin_activity);
                elem_iter->UpdateEdgeMyosinActivity(i,new_edge_myosin_activity);
            }
        }

        // Update SurfaceAreaHistory
        elem_iter->UpdateSurfaceAreaHistory(p_cell_population->rGetMesh().GetSurfaceAreaOfElement(elem_index));

        // Update EdgeLength
        std::vector<double> new_edge_lengths = std::vector<double> (10, 0.6204);
        for (unsigned i = 0; i < elem_iter->GetNumNodes(); ++i)
        {
            unsigned this_node_index = elem_iter->GetNodeGlobalIndex(i);
            unsigned next_node_index = elem_iter->GetNodeGlobalIndex((i + 1) % elem_iter->GetNumNodes());
            double edge_length = p_cell_population->rGetMesh().GetDistanceBetweenNodes(this_node_index, next_node_index);
            new_edge_lengths[i] = edge_length;
        }
        elem_iter->UpdateEdgeLengthHistory(new_edge_lengths);
    }
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyEdgeMyosinActivityModifier)
