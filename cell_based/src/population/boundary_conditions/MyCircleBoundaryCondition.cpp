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

#include "MyCircleBoundaryCondition.hpp"
#include "VertexBasedCellPopulation.hpp"

MyCircleBoundaryCondition::MyCircleBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation,
                                                     c_vector<double, 2> centre,
                                                     double radius)
        : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation),
          mCentre(centre),
          mRadius(radius)
{
    assert(radius>0.0);

    if (dynamic_cast<VertexBasedCellPopulation<2>*>(this->mpCellPopulation)==nullptr)
    {
        EXCEPTION("A VertexBasedCellPopulation must be used with this boundary condition object.");
    }
}

void MyCircleBoundaryCondition::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);
        c_vector<double, 2> node_location = p_node->rGetLocation();
        double signed_distance = norm_2(node_location-mCentre)-mRadius;
        if (signed_distance>0.0)
        {
            c_vector<double, 2> nearest_location = mCentre+mRadius*(node_location-mCentre)/norm_2(node_location-mCentre);
            p_node->rGetModifiableLocation()=nearest_location;
        }
    }
}

bool MyCircleBoundaryCondition::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    unsigned num_nodes=this->mpCellPopulation->GetNumNodes();
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<2>* p_node = this->mpCellPopulation->GetNode(node_index);
        c_vector<double, 2> node_location = p_node->rGetLocation();

        double signed_distance = norm_2(node_location-mCentre)-mRadius;
        if (signed_distance>1e-5)
        {
            condition_satisfied = false;
            break;
        }
    }

    return condition_satisfied;
}

void MyCircleBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyCircleBoundaryCondition)
