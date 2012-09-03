/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "Toroidal2dVertexMesh.hpp"

Toroidal2dVertexMesh::Toroidal2dVertexMesh(double width,
                                           double height,
                                           std::vector<Node<2>*> nodes,
                                           std::vector<VertexElement<2, 2>*> vertexElements,
                                           double cellRearrangementThreshold,
                                           double t2Threshold)
    : MutableVertexMesh<2,2>(nodes, vertexElements, cellRearrangementThreshold, t2Threshold),
      mWidth(width),
      mHeight(height)
{
    // Call ReMesh() to remove any deleted nodes and relabel
    ReMesh();
}

Toroidal2dVertexMesh::Toroidal2dVertexMesh()
{
}

Toroidal2dVertexMesh::~Toroidal2dVertexMesh()
{
}

c_vector<double, 2> Toroidal2dVertexMesh::GetVectorFromAtoB(const c_vector<double, 2>& rLocation1, const c_vector<double, 2>& rLocation2)
{
    assert(mWidth > 0.0);
    assert(mHeight > 0.0);

    c_vector<double, 2> vector = rLocation2 - rLocation1;
    vector[0] = fmod(vector[0], mWidth);
    vector[1] = fmod(vector[1], mHeight);

    // If the points are more than halfway across the domain, measure the other way
    if (vector[0] > 0.5*mWidth)
    {
        vector[0] -= mWidth;
    }
    else if (vector[0] < -0.5*mWidth)
    {
        vector[0] += mWidth;
    }

    // If the points are more than halfway up the domain, measure the other way
    if (vector[1] > 0.5*mHeight)
    {
        vector[1] -= mHeight;
    }
    else if (vector[1] < -0.5*mHeight)
    {
        vector[1] += mHeight;
    }
    return vector;
}

void Toroidal2dVertexMesh::SetNode(unsigned nodeIndex, ChastePoint<2> point)
{
    double x_coord = point.rGetLocation()[0];
    double y_coord = point.rGetLocation()[1];

    // Perform a periodic movement if necessary
    if (x_coord >= mWidth)
    {
        // Move point left
        point.SetCoordinate(0, x_coord - mWidth);
    }
    else if (x_coord < 0.0)
    {
        // Move point right
        point.SetCoordinate(0, x_coord + mWidth);
    }
    if (y_coord >= mHeight)
    {
        // Move point down
        point.SetCoordinate(1, y_coord - mHeight);
    }
    else if (y_coord < 0.0)
    {
        // Move point up
        point.SetCoordinate(1, y_coord + mHeight);
    }

    // Update the node's location
    MutableVertexMesh<2,2>::SetNode(nodeIndex, point);
}

double Toroidal2dVertexMesh::GetWidth(const unsigned& rDimension) const
{
    assert(rDimension==0 || rDimension==1);

    double width = mWidth;
    if (rDimension == 1)
    {
        width = mHeight;
    }

    return width;
}

unsigned Toroidal2dVertexMesh::AddNode(Node<2>* pNewNode)
{
    unsigned node_index = MutableVertexMesh<2,2>::AddNode(pNewNode);

    // If necessary move it to be back onto the torus
    ChastePoint<2> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;
}

double Toroidal2dVertexMesh::GetVolumeOfElement(unsigned index)
{
    VertexElement<2, 2>* p_element = GetElement(index);

    c_vector<double, 2> first_node = p_element->GetNodeLocation(0);
    c_vector<double, 2> current_node;
    c_vector<double, 2> anticlockwise_node;
    c_vector<double, 2> transformed_current_node;
    c_vector<double, 2> transformed_anticlockwise_node;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    double element_area = 0;

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node = p_element->GetNodeLocation(local_index);
        anticlockwise_node = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        /*
         * In order to calculate the area we map the origin to (x[0],y[0])
         * then use GetVectorFromAtoB() to get node coordinates
         */
        transformed_current_node = GetVectorFromAtoB(first_node, current_node);
        transformed_anticlockwise_node = GetVectorFromAtoB(first_node, anticlockwise_node);

        element_area += 0.5*(transformed_current_node[0]*transformed_anticlockwise_node[1]
                             - transformed_anticlockwise_node[0]*transformed_current_node[1]);
    }

    // We take the absolute value just in case the nodes were really oriented clockwise
    return fabs(element_area);
}

c_vector<double, 2> Toroidal2dVertexMesh::GetCentroidOfElement(unsigned index)
{
    VertexElement<2, 2>* p_element = GetElement(index);

    c_vector<double, 2> centroid;
    c_vector<double, 2> transformed_centroid = zero_vector<double>(2);
    c_vector<double, 2> first_node = p_element->GetNodeLocation(0);
    c_vector<double, 2> current_node_location;
    c_vector<double, 2> next_node_location;
    c_vector<double, 2> transformed_current_node;
    c_vector<double, 2> transformed_anticlockwise_node;

    double temp_centroid_x = 0;
    double temp_centroid_y = 0;

    unsigned num_nodes_in_element = p_element->GetNumNodes();

    for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
    {
        // Find locations of current node and anticlockwise node
        current_node_location = p_element->GetNodeLocation(local_index);
        next_node_location = p_element->GetNodeLocation((local_index+1)%num_nodes_in_element);

        /*
         * In order to calculate the centroid we map the origin to (x[0],y[0])
         * then use  GetVectorFromAtoB() to get node coordinates
         */

        transformed_current_node = GetVectorFromAtoB(first_node, current_node_location);
        transformed_anticlockwise_node = GetVectorFromAtoB(first_node, next_node_location);

        temp_centroid_x += (transformed_current_node[0]+transformed_anticlockwise_node[0])*(transformed_current_node[0]*transformed_anticlockwise_node[1]-transformed_current_node[1]*transformed_anticlockwise_node[0]);
        temp_centroid_y += (transformed_current_node[1]+transformed_anticlockwise_node[1])*(transformed_current_node[0]*transformed_anticlockwise_node[1]-transformed_current_node[1]*transformed_anticlockwise_node[0]);
    }

    double vertex_area = GetVolumeOfElement(index);
    double centroid_coefficient = 1.0/(6.0*vertex_area);

    transformed_centroid(0) = centroid_coefficient*temp_centroid_x;
    transformed_centroid(1) = centroid_coefficient*temp_centroid_y;

    centroid = transformed_centroid + first_node;

    return centroid;
}

MutableVertexMesh<2, 2>* Toroidal2dVertexMesh::GetMeshForVtk()
{
    unsigned num_nodes = GetNumNodes();

    std::vector<Node<2>*> temp_nodes(4*num_nodes);
    std::vector<VertexElement<2, 2>*> elements;

    // Create four copies of each node
    for (unsigned index=0; index<num_nodes; index++)
    {
        c_vector<double, 2> location = GetNode(index)->rGetLocation();

        // Node copy at original location
        Node<2>* p_node = new Node<2>(index, false, location[0], location[1]);
        temp_nodes[index] = p_node;

        // Node copy shifted right
        p_node = new Node<2>(num_nodes + index, false, location[0] + mWidth, location[1]);
        temp_nodes[num_nodes + index] = p_node;

        // Node copy shifted up
        p_node = new Node<2>(2*num_nodes + index, false, location[0], location[1] + mHeight);
        temp_nodes[2*num_nodes + index] = p_node;

        // Node copy shifted right and up
        p_node = new Node<2>(3*num_nodes + index, false, location[0] + mWidth, location[1] + mHeight);
        temp_nodes[3*num_nodes + index] = p_node;
    }

    // Iterate over elements
    for (VertexMesh<2,2>::VertexElementIterator elem_iter = GetElementIteratorBegin();
         elem_iter != GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        unsigned num_nodes_in_elem = elem_iter->GetNumNodes();

        std::vector<Node<2>*> elem_nodes;

        // Iterate over nodes contained in this element
        c_vector<double, 2> this_node_location = elem_iter->GetNode(num_nodes_in_elem-1)->rGetLocation();
        for (unsigned local_index=0; local_index<num_nodes_in_elem; local_index++)
        {
            c_vector<double, 2> next_node_location = elem_iter->GetNode(local_index)->rGetLocation();

            unsigned next_node_index = elem_iter->GetNodeGlobalIndex(local_index);

            // Work out whether to use one of the new nodes
            c_vector<double, 2> vector = next_node_location - this_node_location;
            if (vector[0] < -0.5*mWidth)
            {
                next_node_index += num_nodes;
            }
            if (vector[1] < -0.5*mHeight)
            {
                next_node_index += 2*num_nodes;
            }

            elem_nodes.push_back(temp_nodes[next_node_index]);
            this_node_location = temp_nodes[next_node_index]->rGetLocation();
        }

        VertexElement<2,2>* p_element = new VertexElement<2,2>(elem_index, elem_nodes);
        elements.push_back(p_element);
    }

    // Now delete any nodes from the mesh for VTK that are not contained in any elements
    std::vector<Node<2>*> nodes;
    for (unsigned index=0; index<temp_nodes.size(); index++)
    {
        unsigned num_elems_containing_this_node = temp_nodes[index]->rGetContainingElementIndices().size();

        if (num_elems_containing_this_node == 0)
        {
            // Avoid memory leak
            delete temp_nodes[index];
        }
        else
        {
            nodes.push_back(temp_nodes[index]);
        }
    }
    MutableVertexMesh<2, 2>* p_mesh = new MutableVertexMesh<2,2>(nodes, elements, mCellRearrangementThreshold, mT2Threshold);
    return p_mesh;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Toroidal2dVertexMesh)
