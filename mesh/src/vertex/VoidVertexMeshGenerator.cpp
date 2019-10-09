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

#include "VoidVertexMeshGenerator.hpp"

VoidVertexMeshGenerator::VoidVertexMeshGenerator(unsigned numElementsAcross,
                                                           unsigned numElementsUp,
                                                           bool isFlatBottom,
                                                           double cellRearrangementThreshold,
                                                           double t2Threshold,
                                                           double elementArea)
{
    assert(numElementsAcross > 0);
    assert(numElementsUp > 0);
    assert(cellRearrangementThreshold > 0.0);
    assert(t2Threshold > 0.0);
    assert(elementArea > 0.0);

    std::vector<Node<2>*> nodes;
    std::vector<VertexElement<2,2>*>  elements;

    unsigned node_index = 0;
    unsigned node_indices[6];
    unsigned element_index = 0;

    // Each Circle:
    for (unsigned i = numElementsAcross; i<numElementsUp+1; i++)
    {
        bool is_boundary_node = false;
        if (i==numElementsAcross||i==numElementsUp)
            is_boundary_node = true;

        for (unsigned j=1; j<2*i; j++)
        {
            if (j==1)
            {
                double x_coord = i -0.5;
                double y_coord = 1/sqrt(3)/2;
                Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
                nodes.push_back(p_node);
                node_index++;
            }
            else if (j==2)
            {
                double x_coord = i -1.0;
                double y_coord = 1/sqrt(3);
                Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
                nodes.push_back(p_node);
                node_index++;
            }
            else
            {
                double x_coord1 = nodes[node_index-2]->rGetLocation()[0];
                double y_coord1 = nodes[node_index-2]->rGetLocation()[1];

                Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord1-0.5, y_coord1+sqrt(3)/2);
                nodes.push_back(p_node);
                node_index++;
            }
        }
        //rotate:
        for (unsigned j =2*i; j<6*(2*i-1)+1;j++)
        {
            double r = sqrt(pow(nodes[node_index-(2*i-1)]->rGetLocation()[0],2)+pow(nodes[node_index-(2*i-1)]->rGetLocation()[1],2));
            double cos1 = nodes[node_index-(2*i-1)]->rGetLocation()[0]/r;
            double sin1 = nodes[node_index-(2*i-1)]->rGetLocation()[1]/r;
            double x_coord = r*(0.5*cos1-sqrt(3)/2*sin1);
            double y_coord = r*(0.5*sin1+sqrt(3)/2*cos1);
            Node<2>* p_node = new Node<2>(node_index, is_boundary_node, x_coord, y_coord);
            nodes.push_back(p_node);
            node_index++;
        }

    }

    /*
     * Create the elements. The array node_indices contains the
     * global node indices from bottom, going anticlockwise.
     */
    for (unsigned i = numElementsAcross+1; i<numElementsUp+1; i++)
    {
        unsigned n1;
        unsigned n2;
        unsigned n3;

        if (i == numElementsAcross+1)
        {
            n1 = 0;
            n2 = 6*(2*numElementsAcross-1);
            n3 = 6*(-(i-numElementsAcross+1)+(numElementsAcross+i)*(i-numElementsAcross+1));
        }
        else if (i == numElementsAcross+2)
        {
            n1 = 6*(2*numElementsAcross-1);
            n2 = 6*(-(i-1-numElementsAcross+1)+(numElementsAcross+i-1)*(i-numElementsAcross));
            n3 = 6*(-(i-numElementsAcross+1)+(numElementsAcross+i)*(i-numElementsAcross+1));
        }
        else
        {
            n1 = 6*(-(i-2-numElementsAcross+1)+(numElementsAcross+i-2)*(i-1-numElementsAcross));
            n2 = 6*(-(i-1-numElementsAcross+1)+(numElementsAcross+i-1)*(i-numElementsAcross));
            n3 = 6*(-(i-numElementsAcross+1)+(numElementsAcross+i)*(i-numElementsAcross+1));
        }

        for (unsigned j = 1; j<6*(i-1)+1; j++)
        {
            if (j==1)
            {
                node_indices[0] = n3-1;
                node_indices[1] = n3;
                node_indices[2] = n2+1;
                node_indices[3] = n2+2;
                node_indices[4] = n1+1;
                node_indices[5] = n2;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;

            }
            else if (j==2)
            {
                node_indices[0] = n1+1;
                node_indices[1] = n2+2;
                node_indices[2] = n2+3;
                node_indices[3] = n2+4;
                node_indices[4] = n1+3;
                node_indices[5] = n1+2;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }

                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j <i)
            {
                node_indices[0] = elements[element_index-1]->GetNodeGlobalIndex(0)+2;
                node_indices[1] = elements[element_index-1]->GetNodeGlobalIndex(1)+2;
                node_indices[2] = elements[element_index-1]->GetNodeGlobalIndex(2)+2;
                node_indices[3] = elements[element_index-1]->GetNodeGlobalIndex(3)+2;
                node_indices[4] = elements[element_index-1]->GetNodeGlobalIndex(4)+2;
                node_indices[5] = elements[element_index-1]->GetNodeGlobalIndex(5)+2;
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j==i)
            {
                node_indices[0] = n1+(2*(i-1)-1);
                node_indices[1] = n2+(2*((i-1)+1)-1)-1;
                node_indices[2] = n2+(2*((i-1)+1)-1);
                node_indices[3] = n2+(2*((i-1)+1)-1)+1;
                node_indices[4] = n2+(2*((i-1)+1)-1)+2;
                node_indices[5] = n1+(2*(i-1)-1)+1;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;

            }
            else if (j==i+1)
            {
                node_indices[0] = n1+(2*(i-1)-1)+2;
                node_indices[1] = n1+(2*(i-1)-1)+1;
                node_indices[2] = n2+(2*((i-1)+1)-1)+2;
                node_indices[3] = n2+(2*((i-1)+1)-1)+3;
                node_indices[4] = n2+(2*((i-1)+1)-1)+4;
                node_indices[5] = n1+(2*(i-1)-1)+3;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j<2*(i-1)+1)
            {
                node_indices[0] = elements[element_index-1]->GetNodeGlobalIndex(0)+2;
                node_indices[1] = elements[element_index-1]->GetNodeGlobalIndex(1)+2;
                node_indices[2] = elements[element_index-1]->GetNodeGlobalIndex(2)+2;
                node_indices[3] = elements[element_index-1]->GetNodeGlobalIndex(3)+2;
                node_indices[4] = elements[element_index-1]->GetNodeGlobalIndex(4)+2;
                node_indices[5] = elements[element_index-1]->GetNodeGlobalIndex(5)+2;
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j ==2*(i-1)+1)
            {
                node_indices[0] = n1+(2*(i-1)-1)*2;
                node_indices[1] = n2+(2*((i-1)+1)-1)*2-1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*2;
                node_indices[3] = n2+(2*((i-1)+1)-1)*2+1;
                node_indices[4] = n2+(2*((i-1)+1)-1)*2+2;
                node_indices[5] = n1+(2*(i-1)-1)*2+1;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j==2*(i-1)+2)
            {
                node_indices[0] = n1+(2*(i-1)-1)*2+2;
                node_indices[1] = n1+(2*(i-1)-1)*2+1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*2+2;
                node_indices[3] = n2+(2*((i-1)+1)-1)*2+3;
                node_indices[4] = n2+(2*((i-1)+1)-1)*2+4;
                node_indices[5] = n1+(2*(i-1)-1)*2+3;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j<3*(i-1)+1)
            {
                node_indices[0] = elements[element_index-1]->GetNodeGlobalIndex(0)+2;
                node_indices[1] = elements[element_index-1]->GetNodeGlobalIndex(1)+2;
                node_indices[2] = elements[element_index-1]->GetNodeGlobalIndex(2)+2;
                node_indices[3] = elements[element_index-1]->GetNodeGlobalIndex(3)+2;
                node_indices[4] = elements[element_index-1]->GetNodeGlobalIndex(4)+2;
                node_indices[5] = elements[element_index-1]->GetNodeGlobalIndex(5)+2;
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j ==3*(i-1)+1)
            {
                node_indices[0] = n1+(2*(i-1)-1)*3;
                node_indices[1] = n2+(2*((i-1)+1)-1)*3-1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*3;
                node_indices[3] = n2+(2*((i-1)+1)-1)*3+1;
                node_indices[4] = n2+(2*((i-1)+1)-1)*3+2;
                node_indices[5] = n1+(2*(i-1)-1)*3+1;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j==3*(i-1)+2)
            {
                node_indices[0] = n1+(2*(i-1)-1)*3+2;
                node_indices[1] = n1+(2*(i-1)-1)*3+1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*3+2;
                node_indices[3] = n2+(2*((i-1)+1)-1)*3+3;
                node_indices[4] = n2+(2*((i-1)+1)-1)*3+4;
                node_indices[5] = n1+(2*(i-1)-1)*3+3;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j<4*(i-1)+1)
            {
                node_indices[0] = elements[element_index-1]->GetNodeGlobalIndex(0)+2;
                node_indices[1] = elements[element_index-1]->GetNodeGlobalIndex(1)+2;
                node_indices[2] = elements[element_index-1]->GetNodeGlobalIndex(2)+2;
                node_indices[3] = elements[element_index-1]->GetNodeGlobalIndex(3)+2;
                node_indices[4] = elements[element_index-1]->GetNodeGlobalIndex(4)+2;
                node_indices[5] = elements[element_index-1]->GetNodeGlobalIndex(5)+2;
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j ==4*(i-1)+1)
            {
                node_indices[0] = n1+(2*(i-1)-1)*4;
                node_indices[1] = n2+(2*((i-1)+1)-1)*4-1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*4;
                node_indices[3] = n2+(2*((i-1)+1)-1)*4+1;
                node_indices[4] = n2+(2*((i-1)+1)-1)*4+2;
                node_indices[5] = n1+(2*(i-1)-1)*4+1;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j==4*(i-1)+2)
            {
                node_indices[0] = n1+(2*(i-1)-1)*4+2;
                node_indices[1] = n1+(2*(i-1)-1)*4+1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*4+2;
                node_indices[3] = n2+(2*((i-1)+1)-1)*4+3;
                node_indices[4] = n2+(2*((i-1)+1)-1)*4+4;
                node_indices[5] = n1+(2*(i-1)-1)*4+3;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j<5*(i-1)+1)
            {
                node_indices[0] = elements[element_index-1]->GetNodeGlobalIndex(0)+2;
                node_indices[1] = elements[element_index-1]->GetNodeGlobalIndex(1)+2;
                node_indices[2] = elements[element_index-1]->GetNodeGlobalIndex(2)+2;
                node_indices[3] = elements[element_index-1]->GetNodeGlobalIndex(3)+2;
                node_indices[4] = elements[element_index-1]->GetNodeGlobalIndex(4)+2;
                node_indices[5] = elements[element_index-1]->GetNodeGlobalIndex(5)+2;
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j ==5*(i-1)+1)
            {
                node_indices[0] = n1+(2*(i-1)-1)*5;
                node_indices[1] = n2+(2*((i-1)+1)-1)*5-1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*5;
                node_indices[3] = n2+(2*((i-1)+1)-1)*5+1;
                node_indices[4] = n2+(2*((i-1)+1)-1)*5+2;
                node_indices[5] = n1+(2*(i-1)-1)*5+1;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else if (j==5*(i-1)+2)
            {
                node_indices[0] = n1+(2*(i-1)-1)*5+2;
                node_indices[1] = n1+(2*(i-1)-1)*5+1;
                node_indices[2] = n2+(2*((i-1)+1)-1)*5+2;
                node_indices[3] = n2+(2*((i-1)+1)-1)*5+3;
                node_indices[4] = n2+(2*((i-1)+1)-1)*5+4;
                node_indices[5] = n1+(2*(i-1)-1)*5+3;
                for (unsigned jj = 0; jj<6; jj++)
                {
                    node_indices[jj]-= 1;
                }
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            else
            {
                node_indices[0] = elements[element_index-1]->GetNodeGlobalIndex(0)+2;
                node_indices[1] = elements[element_index-1]->GetNodeGlobalIndex(1)+2;
                node_indices[2] = elements[element_index-1]->GetNodeGlobalIndex(2)+2;
                node_indices[3] = elements[element_index-1]->GetNodeGlobalIndex(3)+2;
                node_indices[4] = elements[element_index-1]->GetNodeGlobalIndex(4)+2;
                node_indices[5] = elements[element_index-1]->GetNodeGlobalIndex(5)+2;
                std::vector<Node<2>*> element_nodes;
                for (unsigned k=0; k<6; k++)
                {
                   element_nodes.push_back(nodes[node_indices[k]]);
                }
                VertexElement<2,2>* p_element = new VertexElement<2,2>(element_index, element_nodes);
                elements.push_back(p_element);
                element_index++;
            }
            std::cout<< '\n' << i << ' ' << j << '\n';
            for (unsigned k=0; k<6;k++)
            {
                std::cout<<elements[element_index-1]->GetNodeGlobalIndex(k)<< ' ';
            }


        }
    }

    mpMesh = new MutableVertexMesh<2,2>(nodes, elements, cellRearrangementThreshold, t2Threshold);

    // Scale the mesh so that each element's area takes the value elementArea
    mpMesh->Scale(sqrt(elementArea*2.0/sqrt(3.0)), sqrt(elementArea*2.0/sqrt(3.0)));
}

VoidVertexMeshGenerator::~VoidVertexMeshGenerator()
{
    delete mpMesh;
}

MutableVertexMesh<2,2>* VoidVertexMeshGenerator::GetMesh()
{
    return mpMesh;
}
