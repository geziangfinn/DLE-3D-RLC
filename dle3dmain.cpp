#include <stdio.h>
#include <stdlib.h>
#include "reader.h"
#include "util.h"
#include "topoparser.h"
#include "DME.h"
#include "DLE.h"
#include "netlist.h"
#include "OutputFiles.h"
#include <iostream>

using namespace std;

/**
 * 因为本作业涉及很多专业概念，解决的问题复杂度高
 * 所以请老师在运行/阅读/评审本代码前务必仔细观看 readme 文件
 * 老师辛苦了！
 */

//const char* benchmark_file_path = "D:\\RLC\\DME-3D-RLC\\benchmarks_3D\\testBench"; //todo benchmark.txt文件的绝对路径
//const char* benchmark_file_path = "D:\\RLC\\DME-3D-RLC\\benchmarks_3D\\4test";
int numpins = 0; //存储引脚点的数量
//const char* netlist_file_path = "D:\\RLC\\DME-3D-RLC\\embedding\\netlist.txt"; //todo netlist.txt文件的绝对路径

void print_tree(TreeNode* root,string outputFilePath)
{
    //cout<<"\n*************"<<outputFilePath<<"********\n";
    ofstream preordertree(outputFilePath+"/preOrder.txt");
    ofstream inordertree(outputFilePath+"/inOrder.txt");
    //fstream preorderLayer(outputFilePath+"layerPreOrder.txt")
    
    std::function<void(TreeNode*)> preOrderTraversal = [&](TreeNode* curNode) {
        if(curNode!=nullptr){
            int curId = curNode->ID;
            preordertree<<curId<<" "<<curNode->layer<<endl;        
            preOrderTraversal(curNode->left_child);
            preOrderTraversal(curNode->right_child);
        }
    };

    std::function<void(TreeNode*)> inOrderTraversal = [&](TreeNode* curNode) {
        if(curNode!=nullptr){
            inOrderTraversal(curNode->left_child);
            int curId = curNode->ID;
            inordertree<<curId<<" "<<curNode->layer<<endl;          
            inOrderTraversal(curNode->right_child);
        }
    };

    preOrderTraversal(root);
    inOrderTraversal(root);

}

int main(int argc, char* argv[]) {

    /**
     * benchmark 文件读取解析，获取 2 * numpins 的二维数组，存储所有点的坐标
     */
    //TODO: 调整文件读取结构以适配数据集   ####### Finished #######
    float** sink_set = read(argv[1]);
    printf("The total number of pins: %d\n", numpins);
    for(int i = 0; i< numpins; i++){
        printf("(%f, %f, %d, %d)\n", *(*(sink_set+i)), *(*(sink_set+i)+1), (int)*(*(sink_set+i)+2), (int)*(*(sink_set+i)+3));
    }// 打印点集

    /**
     * 生成二叉树结构
     */
    //TODO: 修改get_nearest以适配3D电路  ###### Finished ######
    TreeNode* tree = topo_generate(sink_set, numpins); //时钟树(二叉树)，返回值是root节点
    TreeNode* RLCtree = topo_generate(sink_set, numpins);

    //TODO: 在生成抽象树之后使用DLE-3D确定merging nodes的layer，应该只需要遍历非sink的节点，确定他们的layer在第一层或第二层
    DLE_3D_multilayer(tree);
    print_tree(tree,argv[2]);
    
    exit(0);
    
    
    DLE_3D(tree);
    DLE_3D(RLCtree);

    //TODO: 在生成网表文件之前就使用DME-3D实现精确位置排布。
    //Step 1
    //对于每两个子节点，用坐标差计算两点连线的斜率，先考虑将merge point放在这条线上移动。
    //Step 2
    //使用Cost函数计算merge point到两个子节点的cost。如果不满足zero skew，那么用while循环进行位置调整。

    //TODO： 如果超过极限而不能平衡？
    DME_3D(tree, false);
    DME_3D(RLCtree, true);

    construct_Snaking_tree(tree);
    construct_Snaking_tree(RLCtree);

    TreeNode* RCsource;
    TreeNode* endLeft;
    TreeNode* endRight;
    RCsource = construct_null_node();
    RCsource->left_child = tree;
    RCsource->layer = tree->layer;
    endLeft = tree->left_child->left_child;

    TreeNode* RLCsource;
    TreeNode* endLeftRLC;
    TreeNode* endRightRLC;
    RLCsource = construct_null_node();
    RLCsource->left_child = RLCtree;
    RLCsource->layer = RLCtree->layer;
    endLeftRLC = RLCtree->left_child->left_child;


    construct_TSV_tree(RCsource);
    outputFiles(RCsource->left_child, false);

    construct_TSV_tree(RLCsource);
    outputFiles(RLCsource->left_child, true);

    float what = exp((-1)*0.2/0.85);
    printf("%f",what);

    /**
     * 生成网表文件
     */
     //generate_netlist(netlist_file_path, tree);

     /**
      * 画出时钟树 (注意看 readme 文件里面对这一部分运行的要求)
      *
      *
      * 需要转到 embedding 程序
      */
    //system("pause");
    return 0;
}