//
// Created by 10403 on 2023/3/10.
//

#ifndef ZST_DME_OUTPUTFILES_H
#define ZST_DME_OUTPUTFILES_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "util.h" //注意文件位置
#include "topoparser.h"
#include "DME.h"

int nodeIndex = 2;
int scriptCount = 0;
int totalWireLength = 0;
int buffernumber = 0;
int TSVnumber = 0;
float distance_constraint = 1500000;

void refreshCount(){
    scriptCount = 0;
};

void outputFiles(TreeNode* root, bool RLC);
void construct_Snaking_tree(TreeNode* node);
bool construct_Split_tree(TreeNode* node, TreeNode* endNode);
void construct_TSV_tree(TreeNode* node);
bool whether_same_layer(TreeNode* node1, TreeNode* node2);
void insert_internal_node(TreeNode* node, TreeNode* childNode);
void insert_extraWL_node(TreeNode* node, TreeNode* childNode);
void inter_split_node(TreeNode* node, TreeNode* childNode);
void insert_buffer_node(TreeNode* node);

void preOrderTraverse(FILE* file, TreeNode* node, bool write);
void preOrderTraverseSinkNode(FILE* file, TreeNode* node, bool write);
void preOrderTraverseWire(FILE* file, TreeNode* node, bool write);
void preOrderTraverseBuffer(FILE* file, TreeNode* node, bool write);

bool whether_same_layer(TreeNode* node1, TreeNode* node2){
    if(!node1 || !node2)
        return true;
    if(node1->layer == node2->layer)
        return true;
    return false;
}

void insert_internal_node(TreeNode* node, TreeNode* childNode){
    TreeNode* TempNode;

    //首先尝试在x方向插入长度为1的一段tsv，标记新节点
    if(fabs(node->x - childNode->x) > 1){
        float new_x;
        if(node->x < childNode->x)
            TempNode = construct_TSV_node(node, node->x+1, node->y);
        else
            TempNode = construct_TSV_node(node, node->x-1, node->y);
    } else if(fabs(node->y - childNode->y) > 1){
        float new_y;
        if(node->y < childNode->y)
            TempNode = construct_TSV_node(node, node->x, node->y+1);
        else
            TempNode = construct_TSV_node(node, node->x, node->y-1);
    } else{
        node->x -= 1;
        TempNode = construct_TSV_node(node, node->x+1, node->y);
    }

 assert(TempNode != NULL);
    if(node->left_child == childNode){
        TempNode->left_child = node->left_child;
        node->left_child = TempNode;
    }

    if(node->right_child == childNode){
        TempNode->left_child = node->right_child;
        node->right_child = TempNode;
    }
}

void insert_extraWL_node(TreeNode* node, TreeNode* childNode){
    if(childNode->extra_WL == 0)
        return;
    //printf("extra: %d, %d\n", node->ID, childNode->ID);
    TreeNode* TempNode;
    if(node->x < childNode->x){
        TempNode = construct_extraWL_node(node, node->x - childNode->extra_WL/2, node->y);
    }else{
        TempNode = construct_extraWL_node(node, node->x + childNode->extra_WL/2, node->y);
    }
    assert(TempNode != NULL);
    if(node->left_child == childNode){
        TempNode->left_child = node->left_child;
        node->left_child = TempNode;
    }

    if(node->right_child == childNode){
        TempNode->left_child = node->right_child;
        node->right_child = TempNode;
    }
}

void inter_split_node(TreeNode* node, TreeNode* childNode){
    float slope = calc_slope(childNode, node);
    float delta_x = distance_constraint/(1+fabs(slope));
    TreeNode* TempNode;
    if(node->x > childNode->x)
        TempNode = construct_buffer_node(childNode, childNode->x + delta_x, childNode->y + delta_x*slope);
    else
        TempNode = construct_buffer_node(childNode, childNode->x - delta_x, childNode->y - delta_x*slope);
    TempNode->C = 35;

    if(node->left_child == childNode){
        TempNode->left_child = node->left_child;
        node->left_child = TempNode;
    }

    if(node->right_child == childNode){
        TempNode->left_child = node->right_child;
        node->right_child = TempNode;
    }
}

void insert_buffer_node(TreeNode* node){
    if(node->needBuffer != 1)
        return;
    //node->inverter1->inverter2->childNode
    buffernumber += 2;

    TreeNode* inverter1 = construct_null_node();
    TreeNode* inverter2 = construct_null_node();

    inverter1->x = node->x;
    inverter1->y = node->y;
    inverter1->layer = node->layer;
    inverter1->needBuffer = 1;

    inverter2->x = node->x;
    inverter2->y = node->y;
    inverter2->layer = node->layer;
    inverter2->needBuffer = 0;

    inverter2->left_child = node->left_child;
    inverter2->right_child = node->right_child;
    inverter1->left_child = inverter2;
    node->left_child = inverter1;
    node->right_child = NULL;
}

void construct_Snaking_tree(TreeNode* node){
    if(node){
        if(node == NULL)
            return;
        //如果已经是根节点，跳过
        if(!node->left_child && !node->right_child){
            return;
        }
        construct_Snaking_tree(node->left_child);
        construct_Snaking_tree(node->right_child);

        //Handle the detouring wiring
        if(node->left_child->extra_WL != 0)
            insert_extraWL_node(node, node->left_child);
        if(node->right_child->extra_WL != 0)
            insert_extraWL_node(node, node->right_child);
    }
}

bool construct_Split_tree(TreeNode* node, TreeNode* endNode){
    bool finalTag = true;
    if(node){
        if(node == endNode)
            return true;
        if(node == NULL)
            return true;
        //如果已经是根节点，跳过
        if(!node->left_child && !node->right_child){
            return true;
        }
        if(construct_Split_tree(node->left_child, endNode) == false)
            finalTag = false;
        if(construct_Split_tree(node->right_child, endNode) == false)
            finalTag = false;
        //Handle the inter split
        if(node->left_child && calc_manhattan_distance(node, node->left_child) > distance_constraint * 2){
            printf("distance: %f\n", calc_manhattan_distance(node, node->left_child));
            inter_split_node(node, node->left_child);
            printf("after split: %f\n", calc_manhattan_distance(node, node->left_child));
            printf("remaining: %f\n\n", calc_manhattan_distance(node->left_child, node->left_child->left_child));
            finalTag = false;
        }
        if(node->right_child && calc_manhattan_distance(node, node->right_child) > distance_constraint * 2){
            inter_split_node(node, node->right_child);
            finalTag = false;
        }
        return finalTag;
        //insert the buffer node
        //insert_buffer_node(node);
    }
}

void longwire_split(TreeNode* node, TreeNode* endNode){
    if(node == endNode)
        return;
    if(node == NULL)
        return;
    //如果已经是根节点，跳过
    if(!node->left_child && !node->right_child){
        return;
    }
//    longwire_split(node->left_child, endNode);
//    longwire_split(node->right_child, endNode);

    if(node->left_child && calc_manhattan_distance(node, node->left_child) > distance_constraint * 2){
        printf("long wire slit\n");
        printf("distance: %f\n", calc_manhattan_distance(node, node->left_child));
        if(calc_manhattan_distance(node, node->right_child) > distance_constraint){
            inter_split_node(node, node->left_child);
            inter_split_node(node, node->right_child);
        }

        printf("after split: %f\n", calc_manhattan_distance(node, node->left_child));
        printf("remaining: %f\n\n", calc_manhattan_distance(node->left_child, node->left_child->left_child));
    }
    if(node->right_child && calc_manhattan_distance(node, node->right_child) > distance_constraint * 2){
        if(calc_manhattan_distance(node, node->left_child) > distance_constraint){
            inter_split_node(node, node->right_child);
            inter_split_node(node, node->left_child);
        }

    }
}

void construct_TSV_tree(TreeNode* node){
    if(node){
        if(node == NULL)
            return;
        //如果已经是根节点，跳过
        if(!node->left_child && !node->right_child){
            return;
        }
        construct_TSV_tree(node->left_child);
        construct_TSV_tree(node->right_child);
        //Handle the inter split
//        if(calc_manhattan_distance(node, node->left_child) > distance_constraint)
//            inter_split_node(node, node->left_child);

        //insert the TSV node
        if(!whether_same_layer(node, node->left_child))
            insert_internal_node(node, node->left_child);
        else if(!whether_same_layer(node, node->right_child))
            insert_internal_node(node, node->right_child);

        //insert the buffer node
        insert_buffer_node(node);
    }
}

void preOrderTraverse(FILE* file, TreeNode* node, bool write){
    if(node == NULL)
        return;
    if(node->left_child == NULL && node->right_child == NULL)
        return;
    if(write){
        fprintf(file, "%d %.0f %.0f\n", nodeIndex, node->x, node->y);
        node->ID = nodeIndex;
        nodeIndex++;
    }
    scriptCount++;
    preOrderTraverse(file, node->left_child, write);
    preOrderTraverse(file, node->right_child, write);
}

void preOrderTraverseSinkNode(FILE* file, TreeNode* node, bool write){
    if(node == NULL)
        return;
    if(node->left_child == NULL && node->right_child == NULL){
        if(write){
            fprintf(file, "%d %d\n", nodeIndex, node->ID);
            node->ID = nodeIndex;
            nodeIndex++;
        }
        scriptCount++;
    }
    preOrderTraverseSinkNode(file, node->left_child, write);
    preOrderTraverseSinkNode(file, node->right_child, write);
}

void preOrderTraverseWire(FILE* file, TreeNode* node, bool write){
    if(node == NULL)
        return;
    if(node->left_child == NULL && node->right_child == NULL)
        return;
    //Only handle the node without buffer
    if(node->needBuffer != 1){
        int delta_WL;
        if(node->left_child != NULL){
            if(write){
                delta_WL = abs(node->x - node->left_child->x) + abs(node->y - node->left_child->y);
                totalWireLength += delta_WL;
                if(node->left_child->TSV == 1){
                    fprintf(file, "%d %d %d\n", node->ID, node->left_child->ID, 1);
                    TSVnumber += 1;
                }

                else
                    fprintf(file, "%d %d %d\n", node->ID, node->left_child->ID, 0);
            }
            scriptCount++;
        }
        if(node->right_child != NULL){
            if(write){
                delta_WL = abs(node->x - node->right_child->x) + abs(node->y - node->right_child->y);
                totalWireLength += delta_WL;
                if(node->right_child->TSV == 1){
                    fprintf(file, "%d %d %d\n", node->ID, node->right_child->ID, 1);
                    TSVnumber += 1;
                }

                else
                    fprintf(file, "%d %d %d\n", node->ID, node->right_child->ID, 0);
            }
            scriptCount++;
        }

    }
    preOrderTraverseWire(file, node->left_child, write);
    preOrderTraverseWire(file, node->right_child, write);
}

void preOrderTraverseBuffer(FILE* file, TreeNode* node, bool write){
    if(node == NULL)
        return;
    if(node->left_child == NULL && node->right_child == NULL)
        return;
    //Only handle the buffer
    if(node->needBuffer == 1){
        assert(node->right_child == NULL && node->left_child != NULL);
        if(write){
            fprintf(file, "%d %d %d\n", node->ID, node->left_child->ID, 0);
        }
        scriptCount++;
    }
    preOrderTraverseBuffer(file, node->left_child, write);
    preOrderTraverseBuffer(file, node->right_child, write);
}

void outputFiles(TreeNode* root, bool RLC){
    FILE* file;
    if(RLC)
        file = fopen("D:\\RLC\\DME-3D-RLC\\RLC_result", "w+");
    else
        file = fopen("D:\\RLC\\DME-3D-RLC\\RC_result", "w+");
    //写入第一行：source Node
    fprintf(file, "%s", "sourcenode 0 0\n");

    int temp;

    //Part2: Num node
    preOrderTraverse(file, root, false);
    fprintf(file, "%s %d\n", "num node", scriptCount+1);

    //用于插入反相器
    fprintf(file, "%s", "1 0 0 \n");
    preOrderTraverse(file, root, true);
    refreshCount();

    //Part3: SinkNode
    preOrderTraverseSinkNode(file, root, false);
    fprintf(file, "%s %d\n", "num sinknode", scriptCount);

    preOrderTraverseSinkNode(file, root, true);
    refreshCount();

    //Part4: wire
    preOrderTraverseWire(file, root, false);
    fprintf(file, "%s %d\n", "num wire", scriptCount+1);

    fprintf(file, "1 2 0\n");
    preOrderTraverseWire(file, root, true);
    refreshCount();

    //Part5: buffer
    preOrderTraverseBuffer(file, root, false);
    fprintf(file, "%s %d\n", "num buffer", scriptCount+1);

    fprintf(file, "0 1 0\n");
    preOrderTraverseBuffer(file, root, true);
    refreshCount();

    printf("TOTALLY %d, TSV: %d, buffer: %d\n", totalWireLength, TSVnumber, buffernumber+1);
    nodeIndex = 2;
    scriptCount = 0;
    totalWireLength = 0;
    TSVnumber = 0;
    buffernumber = 0;
}





#endif //ZST_DME_OUTPUTFILES_H
