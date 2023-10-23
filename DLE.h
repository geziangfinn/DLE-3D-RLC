#ifndef ZST_DME_DLE_H
#define ZST_DME_DLE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <float.h>
#include "util.h" //注意文件位置
#include "topoparser.h"


//在这个部分，用layer=3表示该节点的层级不确定。
void DLE_3D(TreeNode* root);
void DLE_3D_loop(TreeNode* node);
void NearestAssign(TreeNode* node);

void DLE_3D(TreeNode* root){
    DLE_3D_loop(root);
    root->layer = 1;
    NearestAssign(root);
}

void DLE_3D_loop(TreeNode* node){
    if(node){
        //如果已经是根节点，跳过
        if(!node->left_child && !node->right_child){
            return;
        }
        //自下而上递归
        DLE_3D_loop(node->left_child);
        DLE_3D_loop(node->right_child);

        if(node->left_child->layer == 3 || node->right_child->layer == 3)
            node->layer = 3;
        else if(node->left_child->layer == node->right_child->layer)
            node->layer = node->left_child->layer;
        else
            node->layer = 3;
    }
}

void NearestAssign(TreeNode* node){
    //如果已经是根节点，跳过
    if(!node->left_child && !node->right_child)
        return;
    if(node->left_child->layer == 3)
        node->left_child->layer = node->layer;
    if(node->right_child->layer == 3)
        node->right_child->layer = node->layer;

    NearestAssign(node->left_child);
    NearestAssign(node->right_child);
}



#endif //ZST_DME_DLE_H
