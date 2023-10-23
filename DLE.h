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

void DLE_3D_multilayer(TreeNode* root);
void DLE_3D_loop_multilayer(TreeNode* node);
void NearestAssign_multilayer(TreeNode* node);



void DLE_3D_multilayer(TreeNode* root){
    DLE_3D_loop_multilayer(root);
    root->layer = root->el.first;
    assert(root->father==NULL);
    NearestAssign_multilayer(root);
}

void DLE_3D_loop_multilayer(TreeNode *node)
{
    if(node){

        //如果已经是根节点，跳过
        if(!node->left_child && !node->right_child){
            node->el.first=node->layer;
            node->el.second=node->layer;
            return;
        }
        //自下而上递归
        DLE_3D_loop_multilayer(node->left_child);
        DLE_3D_loop_multilayer(node->right_child);

        int l1=max(node->left_child->el.first,node->right_child->el.first);
        int l2=min(node->left_child->el.second,node->right_child->el.second);
        node->el.first=min(l1,l2);
        node->el.second=max(l1,l2);
    }
}

void NearestAssign_multilayer(TreeNode *node)
{
    if(!node->left_child && !node->right_child)
        return;

    if(node->father)
    {   
        assert(node->el.first<=node->el.second);
        if (node->el.first>node->father->layer)
        {
            node->layer=node->el.first;
        }
        else if(node->el.second<node->father->layer)
        {
            node->layer=node->el.second;
        }
        else{
            node->layer=node->father->layer;
        }
    }

    NearestAssign_multilayer(node->left_child);
    NearestAssign_multilayer(node->right_child);
}

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
