#ifndef ZST_DME_TOPOPARSER_H
#define ZST_DME_TOPOPARSER_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "util.h" //注意文件位置
#define MAX 1000
#define c_w 0.0002
#define c_v 15.48
#define r_w 0.0001
#define r_v 0.035
#define c_w_standard c_w/1000000000000000 // 单位F/nm
#define c_v_standard c_v/1000000000000000 // 单位F
#define l_w 0.0000000000000018 // 单位H/nm
#define l_v 0.000000000013827 // 单位H
#define c_constraint 300 // 单位 fF

/**
 * 首先需要说明的是：
 * generate 函数用于生成构建时钟树(二叉树)的路径数组
 * construct 函数用于根据这个路径数组构建时钟树
 * 不能混淆，举个例子：
 * 若我们最终要生成下面这样一棵时钟树：
 *                  R
 *                /   \
 *              S      T
 *            /  \   /  \
 *          A    B C    D
 *  那么 generate 生成的路径数组就是 {(A B S), (C D T), (S T R)}}
 */

typedef struct{
    TreeNode* mergenode;
    TreeNode* leftnode;
    TreeNode* rightnode;
}Triple;

Triple *path[MAX]; //路径数组，用于时钟树(二叉树)的构建，详见 readme 文件“算法步骤”

float calc_dist(TreeNode* , TreeNode* );
float get_nearest(TreeNode *[] , int , TreeNode *[]);
TreeNode* merge_point( TreeNode*  , TreeNode* , int index);
TreeNode* topo_generate(float **sink_set, int n);
void print_path(Triple* [], int );
TreeNode* construct(Triple* [], int);

//void DME_3D(TreeNode* node);
float Downstream_Capacitance(TreeNode* node);
float Elmore_Delay(TreeNode* node, TreeNode* child);
float calc_wire_R(TreeNode* node1, TreeNode* node2);
float calc_wire_C(TreeNode* node1, TreeNode* node2);
//float calc_x_RC(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge);




float Downstream_Capacitance(TreeNode* node){
    if(node){
        //如果已经是根节点，返回当前C值
        if(!node->left_child || !node->right_child){
            return node->C;
        }
        //如果不是，递归计算，加上中间连线的电容值
        float xl = node->left_child -> x;
        float yl = node->left_child -> y;
        float xr = node->right_child -> x;
        float yr = node->right_child -> y;
        int lLayer = node->left_child -> layer;
        int rLayer = node->right_child -> layer;
        float lCap = Downstream_Capacitance(node->left_child);
        float rCap = Downstream_Capacitance(node->right_child);
        float result = (float)c_w*(fabs(xl-xr)+fabs(yl-yr)) + c_v* fabs(lLayer-rLayer) + lCap+rCap;
        return result;
    }
}

float Elmore_Delay(TreeNode* node, TreeNode* child){
    assert(child->father == node);
    float delay = calc_wire_R(node, child) * (0.5 * calc_wire_C(node, child) + child->C);
    return delay;
}



float calc_wire_R(TreeNode* node1, TreeNode* node2){
    float xl = node1 -> x;
    float yl = node1 -> y;
    float xr = node2 -> x;
    float yr = node2 -> y;
    int lLayer = node1 -> layer;
    int rLayer = node2 -> layer;
    float result = (float)r_w*(fabs(xl-xr)+fabs(yl-yr)) + r_v* fabs(lLayer-rLayer);
    return result;
}

float calc_wire_C(TreeNode* node1, TreeNode* node2){
    float xl = node1 -> x;
    float yl = node1 -> y;
    float xr = node2 -> x;
    float yr = node2 -> y;
    int lLayer = node1 -> layer;
    int rLayer = node2 -> layer;
    float result = (float)c_w*(fabs(xl-xr)+fabs(yl-yr)) + c_v* fabs(lLayer-rLayer);
    return result;
}

/**
*计算两点间曼哈顿距离
*@param left 左节点指针
*@param right 右节点指针
*@return 两节点间距离
*/
//使用NN-3D的cost计算方法
float calc_dist(TreeNode* left, TreeNode* right){
    float xl = left -> x;
    float yl = left -> y;
    float xr = right -> x;
    float yr = right -> y;
    int lLayer = left -> layer;
    int rLayer = right -> layer;
    float lCap = left -> C;
    float rCap = right -> C;
    float result = (float)c_w*(fabs(xl-xr)+fabs(yl-yr)) + c_v* fabs(lLayer-rLayer) + 0.1*(lCap+rCap);
    //printf("Coordinate: x %f %f, y %f %f\n", xl, yl, xr, yr);
    //printf("%f\n", (fabs(xl-xr)+fabs(yl-yr)));
    return result;
    //return (float)fabs(xl-xr)+fabs(yl-yr);
}


/**
*计算一组node中最小间距，recurset是递归数组,n是数组长度。我们将结果存进result，输出最小距离
* @param recur_set 一组node
* @param n recur_set的长度
* @param result 这一组中间距最小的两个点
 * 需要指明的是，我们在函数外部就定义了result数组，而且result数组是存储指针(地址)的，
 * 因此当我们计算出最小间距点对并将其保存在result中时
 * 即便跳出这个函数，那个函数外部的result数组仍然保存着最小间距点对
* @return 最小的两点间距离
*/
float get_nearest(TreeNode *recur_set[]  , int n , TreeNode *result[2]){
    int i,j;
    result[0] = recur_set[0];
    result[1] =recur_set[1];

    float minlen = calc_dist(result[0],result[1]);

    //寻找最短Cost的点，并且返回最小cost
    for (i=0;i<n;i++){ // O(N^2)的查找
        for (j=0;j<n;j++){
            if (equals(recur_set[i],recur_set[j])){ //如果碰到一样的两个点(也就是同一个点)就跳过
                continue;
            }
            else if(calc_dist(recur_set[i] , recur_set[j]) < minlen){
                minlen=calc_dist(recur_set[i] , recur_set[j]);
                result[0]=recur_set[i];
                result[1]=recur_set[j];
            }
        }
    }
    return minlen;
}


/**
*计算中点也就是 merging point
* @param left 左节点的指针
* @param right 右节点的指针
* @return 动态内存中开辟的mergepoint的指针
*/
TreeNode* merge_point( TreeNode* left , TreeNode* right, int index){
    TreeNode* mergepoint = (TreeNode *)malloc(sizeof(TreeNode));
    float mergeCapacitance = 0;
    int setLayer = 1;
    float xl = left -> x;
    float yl = left -> y;
    float xr = right -> x;
    float yr = right -> y;
    if(left -> layer == right -> layer)
        setLayer = left -> layer;
    else{
        if(left -> C >= right -> C)
            setLayer = left -> layer;
        else
            setLayer = right -> layer;
    }
    mergeCapacitance = 0;
    mergepoint = convert_to_node( (left -> x + right ->x)/2, (left -> y + right ->y)/2, setLayer, mergeCapacitance, index);
    return mergepoint;

}


/**
*把child和father关联起来形成可访问的二叉树结构
*@param path 存储二叉树生成的路径的数组，即路径数组(见readme文件中算法步骤部分)
*@param n 路径数组的长度
*@return 整个二叉树最底端的根root的指针
**/
TreeNode* construct(Triple* path[], int n){
    int counter = n-1; //先从路径数组的最后一个元素开始，逆序构建
    Triple* temp; //用于迭代
    TreeNode* root= merge(path[n-1] ->leftnode, path[n-1] ->rightnode, path[n-1] ->mergenode); //先合并路径数组中的最后一个 Triple， 即 root 的 Triple
    path[counter]=path[n-1]; //更新
    path[n-1]=NULL; //从路径数组中删除
    counter--;

    while(counter>=0){
        temp=path[counter]; //用于迭代
        assert(temp != NULL); //防止generate中conscounter出现异常
        TreeNode* corres = find(root,temp ->mergenode); //先从已经构建好的那部分时钟树(二叉树)中找到temp节点即corres
        assert(corres != NULL); //可以证明这个corres一定存在( 参考本代码文件最开头的那个例子 )
        insert_left_child(corres,temp -> leftnode);
        insert_right_child(corres,temp -> rightnode);
        path[counter]=NULL;
        counter--;
    }

    return root;
}


/**
*打印路径数组
*@param path 路径数组
*@param n 路径数组的长度
*@return void
**/
void print_path(Triple* path[], int n){
    int i=0;
    float xl,yl,xr,yr,xm,ym;
    int layerLeft, layerRight, layerMerge, IDLeft, IDRight, IDMerge;
    for (i=0;i<n;i++){
        xl= path[i]->leftnode->x;
        yl= path[i]->leftnode->y;
        layerLeft = path[i]->leftnode->layer;
        IDLeft = path[i]->leftnode->ID;
        xr= path[i]->rightnode->x;
        yr= path[i]->rightnode->y;
        layerRight = path[i]->rightnode->layer;
        IDRight = path[i]->rightnode->ID;
        xm= path[i]->mergenode->x;
        ym= path[i]->mergenode->y;
        layerMerge = path[i]->mergenode->layer;
        IDMerge = path[i]->mergenode->ID;
        printf("left point x=%.0f,y=%.0f,layer=%d, ID=%d \nright point x=%.0f,y=%.0f,layer=%d, ID=%d \nmerge point x=%.0f,y=%.0f,layer=%d, ID=%d\n\n",
               xl,yl,layerLeft,IDLeft,xr,yr,layerRight,IDRight,xm,ym,layerMerge,IDMerge);
    }
    putchar('\n');
}


/**
* 生成所有 merging point，并且存储到路径数组里面
*@param sink_set float格式二维数组指针，代表从benchmark文件中读取的坐标点构成的数组
*@param n sink_set数组的行数，也就是点的个数
*@param path 路径数组
*@return 完成二叉树构建之后二叉树的根的指针 root
**/

TreeNode* topo_generate(float **sink_set, int n){
    TreeNode** recur_array = (TreeNode**)malloc(n*sizeof(TreeNode*)); //递归数组；
    for(int i = 0; i < MAX; i++){
        path[i] = (Triple*)malloc(sizeof(Triple*));
    } //初始化路径数组，每个元素为一个Triple的指针
    int conscounter = 0;  //用于路径数组填充的迭代

    //将坐标数组sink_set建立成树节点的集合
    //存入坐标、layer和负载电容
    for (int i=0; i<n; i++){
        recur_array[i] = convert_to_node(*(*(sink_set+i)),*(*(sink_set+i)+1),(int)*(*(sink_set+i)+2),*(*(sink_set+i)+3), i);
    }

    TreeNode* root =NULL;
    TreeNode* merging_point = NULL;
    int arrsize= n; //用于递归数组长度的迭代
    int IDindex = 0;

    while(arrsize > 1){
        TreeNode *result[2]; //用于存储最小间距点对
        //TODO: 修改get_nearest以适配3D电路，不单纯考虑两点之间距离，增加考虑tsv的cost  ###### Finished #######
        float dmin = get_nearest(recur_array, arrsize, result); //最小间距点对的间距
        printf("minimum NN-cost = %f\n", dmin);
        merging_point = merge_point(result[0],result[1], n + IDindex); //计算merging point节点
        IDindex ++;


        // 从递归数组中删除最小间距点对 result[0] 和 result[1]
        // 并将中间节点加入数组
       for(int i=0;i<arrsize;i++){
            if (equals(recur_array[i],result[0]) || equals(recur_array[i],result[1])){
                //要从递归数组里面和删掉一个点的话，不仅要设置成NULL，还得把arrsize-1, 还得把删去之后的空位给补上
                recur_array[i]=recur_array[arrsize-1]; //先把这个空位用最后一个元素补上
                recur_array[arrsize-1]=NULL; //然后把元素删除
                arrsize--; // 递归数组长度减去一
                i=-1; //这里设置为-1, i++之后又变成0，从头开始再找一边，因为result有两个点，得都删除了
            }
        }

        recur_array[arrsize] = merging_point;//新增merging_point
        arrsize++; // 递归数组长度加一

        //建立Triple并填入路径数组
        Triple* triple = (Triple*)malloc(sizeof(Triple));
        triple -> leftnode = result[0];
        triple -> rightnode = result[1];
        triple ->mergenode = merging_point;

        path[conscounter] = triple;
        conscounter++;

    }
    //检查路径数组是否有问题
    //print_path(path,conscounter);

    //完成二叉树构建
    root = construct(path,conscounter);
    //收回空间
    free(recur_array);

    return root;
}
#endif //ZST_DME_TOPOPARSER_H
