/*

Copyright (c) 2005-2019, University of Oxford.
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
#include "VertexElement.hpp"
#include <cassert>
#include <iostream>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes),
      mFaces(rFaces),
      mOrientations(rOrientations)
{
    // This constructor should only be used in 3D
    assert(SPACE_DIM == 3);    // LCOV_EXCL_LINE - code will be removed at compile time

    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    if (SPACE_DIM == ELEMENT_DIM)
    {
        // Register element with nodes
        this->RegisterWithNodes();
    }
    
    // My changes.
    mSurfaceAreaHistory = std::vector<double> (50, 0.6204*6.0);//0.6204 = Ra

    mEdgeLengthHistory = std::vector<std::vector<double>> (50, std::vector<double> (10, 0.6204));
//    mEdgeMyosinActivities = std::vector<double> (6, 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<VertexElement<ELEMENT_DIM-1,SPACE_DIM>*>& rFaces,
                                                     const std::vector<bool>& rOrientations)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index),
      mFaces(rFaces),
      mOrientations(rOrientations)
{
    // Each face must have an associated orientation
    assert(mFaces.size() == mOrientations.size());

    // Make a set of nodes with mFaces
    std::set<Node<SPACE_DIM>* > nodes_set;
    for (unsigned face_index=0; face_index<mFaces.size(); face_index++)
    {
        for (unsigned node_index=0; node_index<mFaces[face_index]->GetNumNodes(); node_index++)
        {
            nodes_set.insert(mFaces[face_index]->GetNode(node_index));
        }
    }

    // Populate mNodes
    for (typename std::set< Node<SPACE_DIM>* >::iterator node_iter = nodes_set.begin();
         node_iter != nodes_set.end();
         ++node_iter)
    {
         this->mNodes.push_back(*node_iter);
    }

    // Register element with nodes
    this->RegisterWithNodes();
    
    // My changes.
    mSurfaceAreaHistory = std::vector<double> (50, 0.6204*6.0);

    mEdgeLengthHistory = std::vector<std::vector<double>> (50, std::vector<double> (10, 0.6204));
//    mEdgeMyosinActivities = std::vector<double> (6, 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index)
{
    // My changes.
    mSurfaceAreaHistory = std::vector<double> (50, 0.6204*6.0);

    mEdgeLengthHistory = std::vector<std::vector<double>> (50, std::vector<double> (10, 0.6204));
//    mEdgeMyosinActivities = std::vector<double> (6, 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::VertexElement(unsigned index,
                                                     const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<ELEMENT_DIM, SPACE_DIM>(index, rNodes)
{
    // My changes.
    mSurfaceAreaHistory = std::vector<double> (50, 0.6204*6.0);

    mEdgeLengthHistory = std::vector<std::vector<double>> (50, std::vector<double> (10, 0.6204));
//    mEdgeMyosinActivities = std::vector<double> (6, 1.0);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM, SPACE_DIM>::~VertexElement()
{
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned VertexElement<ELEMENT_DIM, SPACE_DIM>::GetNumFaces() const
{
    return mFaces.size();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::AddFace(VertexElement<ELEMENT_DIM-1,SPACE_DIM>* pFace)
{
    // Add pFace to the end of mFaces
    this->mFaces.push_back(pFace);

    // Create a set of indices of nodes currently owned by this element
    std::set<unsigned> node_indices;
    for (unsigned local_index=0; local_index<this->GetNumNodes(); local_index++)
    {
        node_indices.insert(this->GetNodeGlobalIndex(local_index));
    }

    // Loop over nodes owned by pFace
    unsigned end_index = this->GetNumNodes()-1;
    for (unsigned local_index=0; local_index<pFace->GetNumNodes(); local_index++)
    {
        // If this node is not already owned by this element...
        unsigned global_index = pFace->GetNodeGlobalIndex(local_index);
        if (node_indices.find(global_index) == node_indices.end())
        {
            // ... then add it to the element (and vice versa)
            this->AddNode(pFace->GetNode(local_index), end_index);
            end_index++;
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
VertexElement<ELEMENT_DIM-1,  SPACE_DIM>* VertexElement<ELEMENT_DIM, SPACE_DIM>::GetFace(unsigned index) const
{
    assert(index < mFaces.size());
    return mFaces[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool VertexElement<ELEMENT_DIM, SPACE_DIM>::FaceIsOrientatedClockwise(unsigned index) const
{
    assert(index < mOrientations.size());
    return mOrientations[index];
}

// My changes

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetHistoricSurfaceArea() const
{
    return mSurfaceAreaHistory.front();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateSurfaceAreaHistory(double currentSurfaceArea)
{
    mSurfaceAreaHistory.erase(mSurfaceAreaHistory.begin());
    mSurfaceAreaHistory.push_back(currentSurfaceArea);
}

// My changes
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetEdgeMyosinActivity(unsigned index) const
{
//std::vector<double>::size_type index_edge = index;
//if (SimulationTime::Instance()->GetTime() == 31.0)
//std::cout << mEdgeMyosinActivities.size();
//std::cout << index_edge;
    assert(!this-> IsDeleted());
    //assert(index_edge < mEdgeMyosinActivities.size());
    return mEdgeMyosinActivities[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::InitializeEdgeMyosinActivities(double newMyosinActivity)
{
    //mEdgeMyosinActivities = std::vector<double> (this->GetNumNodes(), newMyosinActivity);
    mEdgeMyosinActivities.reserve(this->GetNumNodes());
    for (unsigned i = 0; i < this->GetNumNodes(); i++)
    {
        c_vector<double, SPACE_DIM> point1= this->mNodes[i]->rGetLocation();
        double xcoord1 = point1[0];
        c_vector<double, SPACE_DIM> point2= this->mNodes[(i+1)%this->GetNumNodes()]->rGetLocation();
        double xcoord2 = point2[0];
        double length = norm_2(point2-point1);
        double theta = acos((xcoord2-xcoord1)/length);
        if (theta > M_PI/2)
            theta = M_PI -theta;
        mEdgeMyosinActivities[i] = newMyosinActivity*(0.5 + 2.0/M_PI*theta);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateEdgeMyosinActivity(unsigned index,double newEdgeMyosinActivity)
{
    mEdgeMyosinActivities[index] = newEdgeMyosinActivity;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::InitializeEdgeLengthHistory()
{
    mEdgeLengthHistory = std::vector<std::vector<double>> (50, std::vector<double> (10, 0.6204));
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void VertexElement<ELEMENT_DIM, SPACE_DIM>::UpdateEdgeLengthHistory(std::vector<double> newEdgeLengths)
{
    mEdgeLengthHistory.erase(mEdgeLengthHistory.begin());
    mEdgeLengthHistory.push_back(newEdgeLengths);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetHistoricEdgeLength(unsigned index) const
{
    //std::vector<double> historic_edge_lengths = mEdgeLengthHistory.front();
    //return historic_edge_lengths[index];
    return mEdgeLengthHistory.front()[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double VertexElement<ELEMENT_DIM, SPACE_DIM>::GetEdgeLength(unsigned index) const
{
    std::vector<double> vec = *(--mEdgeLengthHistory.end());
    return vec[index];
}

//////////////////////////////////////////////////////////////////////
//                  Specialization for 1d elements                  //
//                                                                  //
//                 1d elements are just edges (lines)               //
//////////////////////////////////////////////////////////////////////

/**
 * Specialization for 1d elements so we don't get errors from Boost on some
 * compilers.
 */
template<unsigned SPACE_DIM>
VertexElement<1, SPACE_DIM>::VertexElement(unsigned index, const std::vector<Node<SPACE_DIM>*>& rNodes)
    : MutableElement<1, SPACE_DIM>(index, rNodes)
{
}

template<unsigned SPACE_DIM>
unsigned VertexElement<1, SPACE_DIM>::GetNumFaces() const
{
    return 0;
}

template<unsigned SPACE_DIM>
VertexElement<0, SPACE_DIM>* VertexElement<1, SPACE_DIM>::GetFace(unsigned index) const
{
    return nullptr;
}

template<unsigned SPACE_DIM>
bool VertexElement<1, SPACE_DIM>::FaceIsOrientatedClockwise(unsigned index) const
{
    return false;
}

// My changes
template<unsigned SPACE_DIM>
double VertexElement<1, SPACE_DIM>::GetHistoricSurfaceArea() const
{
    return 0.0;
}

template<unsigned SPACE_DIM>
double VertexElement<1, SPACE_DIM>::GetEdgeMyosinActivity(unsigned index) const
{
    return 0.0;
}

template<unsigned SPACE_DIM>
double VertexElement<1, SPACE_DIM>::GetHistoricEdgeLength(unsigned index) const
{
    return 0.0;
}

template<unsigned SPACE_DIM>
double VertexElement<1, SPACE_DIM>::GetEdgeLength(unsigned index) const
{
    return 0.0;
}

// Explicit instantiation
template class VertexElement<1,1>;
template class VertexElement<1,2>;
template class VertexElement<1,3>;
template class VertexElement<2,2>;
template class VertexElement<2,3>;
template class VertexElement<3,3>;
