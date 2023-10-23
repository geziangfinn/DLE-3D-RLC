//
// Created by 10403 on 2023/2/16.
//

#ifndef ZST_DME_DME_H
#define ZST_DME_DME_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <stdbool.h>
#include <float.h>
#include "util.h" //注意文件位置
#include "topoparser.h"
float skewBoundary = 0.01;
float skewModifyStep = 1;

struct coordinate{
    float x, y;
};

typedef struct coordinate coordinate;

bool cmp(coordinate a,coordinate b)
{
    if(a.x==b.x)
        return a.y<b.y;
    return a.x<b.x;
}



bool whether_leaf(TreeNode* node);
float calc_manhattan_distance(TreeNode* a, TreeNode* b);
float calc_standard_Capacitance(float capacitance_in_fF);
float calc_slope(TreeNode* a, TreeNode* b);
void update_merge_coord(TreeNode* merge, TreeNode* a, float slope, float manhattan_distance);
void update_merge_Capacitance(TreeNode* nodeMerge, TreeNode * nodeLeft, TreeNode* nodeRight, float ea, float eb);
void update_merge_Delay(TreeNode* nodeMerge, TreeNode * nodeLeft, TreeNode* nodeRight, float ea, float eb);
void update_merge_Delay_RLC(TreeNode* nodeMerge, TreeNode * nodeLeft, TreeNode* nodeRight, float ea, float eb);
float calc_delay_RLC(TreeNode* nodeMerge, TreeNode* nodeChild, float WL);
float calc_delay_RC(TreeNode* nodeMerge, TreeNode* nodeChild, float WL);
void modify_coord_by_L1(float* ea_pointer, float* eb_pointer, TreeNode* nodeLeft, TreeNode* nodeRight, float x);//For the situation that delay_a > delay_b in RLC
void modify_coord_by_L2(float* ea_pointer, float* eb_pointer, TreeNode* nodeLeft, TreeNode* nodeRight, float x);//For the situation that delay_a < delay_b in RLC

void DME_3D(TreeNode* node, bool RLC);
void calculation(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge, bool RLC);
//在RC模型下计算步骤中的线长x
float calc_x_RC(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge, float L);
float calc_L2_RC(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge, int tag);

void RC_calculation(TreeNode* nodeMerge, TreeNode* nodeLeft, TreeNode* nodeRight, float* ea_pointer, float* eb_pointer, float L, float slope);
void RLC_calculation(TreeNode* nodeMerge, TreeNode* nodeLeft, TreeNode* nodeRight, float* ea_pointer, float* eb_pointer, float L, float slope);

bool whether_leaf(TreeNode* node){
    if(node->left_child == NULL && node->right_child == NULL)
        return true;
    return false;
}

float calc_manhattan_distance(TreeNode* a, TreeNode* b){
    return (abs(a->x - b->x) + abs(a->y - b->y));
}

float calc_standard_Capacitance(float capacitance_in_fF){
    return (capacitance_in_fF/1000000000000000);
}

float calc_slope(TreeNode* a, TreeNode* b){
    return ((b->y-a->y) / (b->x-a->x));
}

void update_merge_coord(TreeNode* merge, TreeNode* a, float slope, float manhattan_distance){
    float delta_x, delta_y;
    delta_x = manhattan_distance * (1 / (1+fabs(slope) ) );
    delta_y = delta_x * slope;
    merge->x = a->x + delta_x;
    merge->y = a->y + delta_y;
}

void update_merge_Capacitance(TreeNode* nodeMerge, TreeNode * nodeLeft, TreeNode* nodeRight, float ea, float eb){
    float delta_C = nodeLeft->C + nodeRight->C + c_w*(ea + eb) + c_v*(abs(nodeRight->layer - nodeLeft->layer));
    //考虑了buffer insertion的电容update
    if(delta_C > c_constraint){
        nodeMerge->C = 300;
        nodeMerge->needBuffer = 1;
        return;
    }
    nodeMerge->C = delta_C;
}

void update_merge_Delay(TreeNode* nodeMerge, TreeNode * nodeLeft, TreeNode* nodeRight, float ea, float eb){
    float delta_delay = nodeLeft->delay + 0.695 * (
            0.5 * r_w * c_w * ea * ea +
            r_w * (nodeLeft->C + c_v * abs(nodeMerge->layer - nodeLeft->layer)) * ea +
            r_v * (nodeLeft->C * abs(nodeMerge->layer - nodeLeft->layer) + 0.5 * c_v * (nodeMerge->layer - nodeLeft->layer) * (nodeMerge->layer - nodeLeft->layer)) -
            (r_w * c_v - r_v * c_w) * abs(nodeMerge->layer - nodeLeft->layer) * ea);
    float delta_delay_right = nodeRight->delay + 0.695 * (
                              0.5 * r_w * c_w * eb * eb +
                              r_w * (nodeRight->C + c_v * abs(nodeMerge->layer - nodeRight->layer)) * eb +
                              r_v * (nodeRight->C * abs(nodeMerge->layer - nodeRight->layer) + 0.5 * c_v * (nodeMerge->layer - nodeRight->layer) * (nodeMerge->layer - nodeRight->layer)) -
                              (r_w * c_v - r_v * c_w) * abs(nodeMerge->layer - nodeRight->layer) * eb);

    /*float delta_delay_right = nodeRight->delay +
            r_w * (ea + eb) * (nodeRight->C + 0.5 * c_w * (ea + eb)) +
            0.5 * r_w * c_w * ea * ea -
            r_w * (nodeRight->C + c_v * abs(nodeMerge->layer - nodeRight->layer) + c_w * (ea + eb)) * ea +
            r_v * (nodeRight->C * abs(nodeMerge->layer - nodeRight->layer) + 0.5 * c_v * abs(nodeMerge->layer - nodeRight->layer) * abs(nodeMerge->layer - nodeRight->layer)) +
            r_w * c_w * abs(nodeMerge->layer - nodeRight->layer) * (ea + eb) -
            (r_w * c_v - r_v * c_w) * abs(nodeMerge->layer - nodeRight->layer) * ea;*/
    printf("/****Test Merge Delay****/\n"
           "Delay before merge: Left: %f, Right: %f, after: Left: %f,  Right: %f\n", nodeLeft->delay,nodeRight->delay, delta_delay, delta_delay_right);
    //按最大delay进行录入
    if(delta_delay >= delta_delay_right)
        nodeMerge->delay = delta_delay;
    else
        nodeMerge->delay = delta_delay_right;
}

void update_merge_Delay_RLC(TreeNode* nodeMerge, TreeNode * nodeLeft, TreeNode* nodeRight, float ea, float eb){
    float delta_delay_left, delta_delay_right;
    float var_left, var_right;
    var_left = calc_delay_RLC(nodeMerge, nodeLeft, ea);
    var_right = calc_delay_RLC(nodeMerge, nodeRight, eb);

    delta_delay_left = nodeLeft->delay + var_left;
    delta_delay_right = nodeRight->delay + var_right;

    printf("/****Test RLC Merge Delay****/\n"
           "Delay before merge: Left: %f, Right: %f, after: Left: %f,  Right: %f\n", nodeLeft->delay,nodeRight->delay, delta_delay_left, delta_delay_right);
    //按最大delay进行录入
    if(delta_delay_left >= delta_delay_right)
        nodeMerge->delay = delta_delay_left;
    else
        nodeMerge->delay = delta_delay_right;
}

float calc_x_RC(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge, float L){
    float beta = r_v * (abs(nodeRight->layer - nodeMerge->layer)*nodeRight->C -
                        abs(nodeMerge->layer - nodeLeft->layer)*nodeLeft->C +
                        0.5*c_v*((nodeRight->layer - nodeMerge->layer)^2 -
                        (nodeMerge->layer - nodeLeft->layer)^2));
    float up = (nodeRight->delay - nodeLeft->delay) + r_w * L * (nodeRight->C + 0.5 * c_w * L) + beta + r_v * c_w * abs(nodeRight->layer - nodeMerge->layer) * L;
    float down = r_w * (nodeLeft->C + nodeRight->C + c_w * L) + r_v * c_w * abs(nodeRight->layer - nodeLeft->layer);
    float x = up/down;
    return x;
}

float calc_L2_RC(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge, int tag){
    float alpha, beta, up;
    //tag = 0: |eb| = L'
    //tag = 1: |ea| = L'
    if(tag == 0){
        alpha = r_w * nodeRight->C + r_v * c_w * abs(nodeRight->layer - nodeMerge->layer);
        beta = r_v * (abs(nodeRight->layer - nodeMerge->layer)*nodeRight->C -
                      abs(nodeMerge->layer - nodeLeft->layer)*nodeLeft->C +
               0.5*c_v*((nodeRight->layer - nodeMerge->layer)^2 -
                                                              (nodeMerge->layer - nodeLeft->layer)^2));
        up = sqrt(2 * r_w * c_w * (nodeLeft->delay - nodeRight->delay) + alpha * alpha - 2 * r_w * c_w * beta) - alpha;
    }

    else{
        alpha = r_w * nodeLeft->C + r_v * c_w * abs(nodeLeft->layer - nodeMerge->layer);
        beta = r_v * (abs(nodeMerge->layer - nodeLeft->layer)*nodeLeft->C -
                      abs(nodeRight->layer - nodeMerge->layer)*nodeRight->C +
               0.5*c_v*((nodeRight->layer - nodeMerge->layer)^2 -
                                                              (nodeMerge->layer - nodeLeft->layer)^2));
        up = sqrt(2 * r_w * c_w * (nodeRight->delay - nodeLeft->delay) + alpha * alpha - 2 * r_w * c_w * beta) - alpha;
    }
    return up/(r_w * c_w);
}

void RC_calculation(TreeNode* nodeMerge, TreeNode* nodeLeft, TreeNode* nodeRight, float* ea_pointer, float* eb_pointer, float L, float slope){
    float t_a, t_b;
    t_a = calc_delay_RC(nodeMerge, nodeLeft, *ea_pointer);
    t_b = calc_delay_RC(nodeMerge, nodeRight, *eb_pointer);
    while(fabs(t_a + nodeLeft->delay - t_b - nodeRight->delay) >= 1){
        if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
            modify_coord_by_L1(ea_pointer, eb_pointer, nodeLeft, nodeRight, skewModifyStep);
        else
            modify_coord_by_L2(ea_pointer, eb_pointer, nodeLeft, nodeRight, skewModifyStep);

        t_a = calc_delay_RC(nodeMerge, nodeLeft, *ea_pointer);
        t_b = calc_delay_RC(nodeMerge, nodeRight, *eb_pointer);
    }
    printf("left: %f, right: %f\n t_a: %f, t_b: %f\n total left: %f, total right: %f\n", nodeLeft->delay, nodeRight->delay, t_a, t_b, t_a+nodeLeft->delay, t_b + nodeRight->delay);
    printf("RC delay: left: %f, right: %f\n", t_a, t_b);
    if(nodeLeft->extra_WL!=0)
        nodeLeft->extra_WL = *ea_pointer - L;
    else if(nodeRight->extra_WL!=0)
        nodeRight->extra_WL = *eb_pointer - L;
    else
        update_merge_coord(nodeMerge, nodeLeft, slope, *ea_pointer);
    update_merge_Capacitance(nodeMerge, nodeLeft, nodeRight, *ea_pointer, *eb_pointer);
    //update_merge_Delay(nodeMerge, nodeLeft, nodeRight, ea, eb);
    if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
        nodeMerge->delay = t_a + nodeLeft->delay;
    else
        nodeMerge->delay = t_b + nodeRight->delay;
}


void RLC_calculation(TreeNode* nodeMerge, TreeNode* nodeLeft, TreeNode* nodeRight, float* ea_pointer, float* eb_pointer, float L, float slope){
    float t_a, t_b;
    t_a = calc_delay_RLC(nodeMerge, nodeLeft, *ea_pointer);
    t_b = calc_delay_RLC(nodeMerge, nodeRight, *eb_pointer);
    //需要考虑extra wirelength存在的情况。但是尽量不要引入额外的线长。
    while(fabs(t_a + nodeLeft->delay - t_b - nodeRight->delay) >= 1){
        if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
            modify_coord_by_L1(ea_pointer, eb_pointer, nodeLeft, nodeRight, skewModifyStep);
        else
            modify_coord_by_L2(ea_pointer, eb_pointer, nodeLeft, nodeRight, skewModifyStep);

        t_a = calc_delay_RLC(nodeMerge, nodeLeft, *ea_pointer);
        t_b = calc_delay_RLC(nodeMerge, nodeRight, *eb_pointer);
    }
    //printf("原本的 delay: left: %f, right: %f\n", nodeLeft->delay, nodeRight->delay);
    printf("left: %f, right: %f\n t_a: %f, t_b: %f\n total left: %f, total right: %f\n", nodeLeft->delay, nodeRight->delay, t_a, t_b, t_a+nodeLeft->delay, t_b + nodeRight->delay);
    t_a = calc_delay_RC(nodeMerge, nodeLeft, *ea_pointer);
    t_b = calc_delay_RC(nodeMerge, nodeRight, *eb_pointer);
    printf("RC delay: %f, %f\n", t_a, t_b);
    if(nodeLeft->extra_WL!=0 && *ea_pointer > L)
        nodeLeft->extra_WL = *ea_pointer - L;
    else if(nodeRight->extra_WL!=0 && *eb_pointer > L)
        nodeRight->extra_WL = *eb_pointer - L;
    else
        update_merge_coord(nodeMerge, nodeLeft, slope, *ea_pointer);
    update_merge_Capacitance(nodeMerge, nodeLeft, nodeRight, *ea_pointer, *eb_pointer);
    //update_merge_Delay(nodeMerge, nodeLeft, nodeRight, ea, eb);
    if(t_a + nodeLeft->delay > t_b + nodeRight->delay)
        nodeMerge->delay = t_a + nodeLeft->delay;
    else
        nodeMerge->delay = t_b + nodeRight->delay;
}

float calc_delay_RLC(TreeNode* nodeMerge, TreeNode* nodeChild, float WL){
    float t_pdi, theta, omega, numerator, denominator, elmore;
    float wireLength = WL + nodeChild->extra_WL;
    if(wireLength == 0)
        return 0;
    //如果两节点中间没有TSV
    assert(nodeChild->C <= c_constraint);
    numerator = wireLength * r_w * (0.5 * wireLength * c_w_standard + calc_standard_Capacitance(nodeChild->C));
    //这里分母还没算完，后续要开根
    denominator = wireLength * l_w * (0.5 * wireLength * c_w_standard + calc_standard_Capacitance(nodeChild->C));
    if(nodeMerge->layer != nodeChild->layer){
        numerator += r_v * (0.5 * c_v_standard + calc_standard_Capacitance(nodeChild->C + wireLength * c_w));
        //denominator += l_v * (0.5 * c_v_standard + calc_standard_Capacitance(nodeChild->C + wireLength * c_w));
    }
    //给分母开根
    denominator = sqrt(denominator);
    //154953990144.000000
    theta = 0.5 * (numerator/denominator);
    omega = 1/denominator;
    elmore = 0.695 * numerator;
    //将单位换算回ps
    t_pdi = roundf(1000000000000000 * (  (1.047 * exp((-1)*theta/0.85))/omega + elmore  ));
    //printf("before : %f\n", 1000000000000000 * (  (1.047 * exp((-1)*theta/0.85))/omega));
    return t_pdi;
}

//这里算的是到该节点的delay，不应该用这个作为指标进行delay平衡！！
float calc_delay_RC(TreeNode* nodeMerge, TreeNode* nodeChild, float WL){
    float t_pdi, theta, omega, numerator, denominator, elmore;
    float wireLength = WL + nodeChild->extra_WL;
    //如果两节点中间没有TSV
    numerator = wireLength * r_w * (0.5 * wireLength * c_w_standard + calc_standard_Capacitance(nodeChild->C));
    if(nodeMerge->layer != nodeChild->layer)
        numerator += r_v * (0.5 * c_v_standard + calc_standard_Capacitance(nodeChild->C + wireLength * c_w));
    elmore = 0.695 * numerator;
    t_pdi = 1000000000000000 * elmore;
    return t_pdi;
}

//delay_a > delay_b, ea -= x/2, eb+= x/2
void modify_coord_by_L1(float* ea_pointer, float* eb_pointer, TreeNode* nodeLeft, TreeNode* nodeRight, float x){
    float delta_x;
    //if it will cause nodeRight to have extra wirelength
    if(*ea_pointer < x){
        delta_x = x - *ea_pointer;
        *ea_pointer = 0;
        *eb_pointer += delta_x;
        nodeRight->extra_WL += delta_x;
        return;
    }
    //if nodeLeft already have extra wirelength
    else if(nodeLeft->extra_WL > 0){
        if(nodeLeft->extra_WL < x){
            delta_x = x - nodeLeft->extra_WL;
            nodeLeft->extra_WL = 0;
            *ea_pointer -= delta_x/2;
            *eb_pointer += delta_x/2;
        } else{
            nodeLeft->extra_WL -= x;
            *ea_pointer -= x;
        }
    }
    //if ea still in [0, L] range
    else{
        *ea_pointer -= x/2;
        *eb_pointer += x/2;
    }
}

//delay_a < delay_b, ea += x/2, eb -= x/2
void modify_coord_by_L2(float* ea_pointer, float* eb_pointer, TreeNode* nodeLeft, TreeNode* nodeRight, float x){
    float delta_x;
    //if it will cause nodeLeft to have extra wirelength
    if(*eb_pointer < x){
        delta_x = x - *eb_pointer;
        *eb_pointer = 0;
        *ea_pointer += delta_x;
        nodeLeft->extra_WL += delta_x;
        return;
    }
        //if nodeRight already have extra wirelength
    else if(nodeRight->extra_WL > 0){
        if(nodeRight->extra_WL < x){
            delta_x = x - nodeRight->extra_WL;
            nodeRight->extra_WL = 0;
            *ea_pointer += delta_x/2;
            *eb_pointer -= delta_x/2;
        } else{
            nodeRight->extra_WL -= x;
            *eb_pointer -= x;
        }
    }
        //if eb still in [0, L] range
    else{
        *ea_pointer += x/2;
        *eb_pointer -= x/2;
    }
}

void calculation(TreeNode* nodeLeft, TreeNode* nodeRight, TreeNode* nodeMerge, bool RLC){
    //该函数假定nodeLeft的x小于nodeRight的x
    assert(nodeLeft->x <= nodeRight->x);

    //先假设能在0<=x<=k的范围内有解。
    //第一步 求出d(ms(a), ms(b))=k, 在这里k=L。单独使用函数进行判断,已经在DME.h中实现。
    float L = calc_manhattan_distance(nodeLeft, nodeRight);

    //第二步，嵌套公式求出实际的merging cost，并与L进行比较
    float x = calc_x_RC(nodeLeft, nodeRight, nodeMerge, L);

    //第三步，比较x是否在0到L的区间内，由此得到ea,eb
    float ea, eb, slope;
    slope = 0;
    //如果算出来的x在0到L的范围内，那么有解，|ea| = x, |eb| = L - x
    //找斜率，确定merge node的精确坐标
    //assert(0<= x && x <= L);
    if(0 <= x && x <= L){
        ea = x;
        eb = L - x;
        printf("ea: %f, eb: %f ", ea, eb);
        slope = calc_slope(nodeLeft, nodeRight);
        printf("slope: %f\n", slope);
        update_merge_coord(nodeMerge, nodeLeft, slope, ea);
    }
    //如果小于0：那么|ea| = 0, |eb| = L', 反之则是|ea| = L', |eb| = 0
    //这个情况的意思应该是，merge node与其中一个重合，对于另一个节点则需要做snaking
    else if(x < 0){
        printf("now<0\n");
        ea = 0;
        eb = calc_L2_RC(nodeLeft, nodeRight, nodeMerge, 0);
        assert(eb > 0);
        //更新子节点的extra_WL，更新merge node的坐标
        nodeRight->extra_WL = eb - L;
        nodeMerge->x = nodeLeft->x;
        nodeMerge->y = nodeLeft->y;
    }
    else if(x > L){
        printf("now>L\n");
        printf("left: %d, right: %d\n", nodeLeft->ID, nodeRight->ID);
        ea = calc_L2_RC(nodeLeft, nodeRight, nodeMerge, 1);
        assert(ea > 0);
        eb = 0;
        nodeLeft->extra_WL = ea - L;
        nodeMerge->x = nodeRight->x;
        nodeMerge->y = nodeRight->y;
    }
    nodeMerge->x = round(nodeMerge->x);
    nodeMerge->y = round(nodeMerge->y);
    //第四步，更新Merge node的C
    //即便引入extra_WL,该函数也不需改动
    update_merge_Capacitance(nodeMerge, nodeLeft, nodeRight, ea, eb);

    //第五步，更新Merge node的delay
    //即便引入extra_WL,该函数也不需改动
    update_merge_Delay(nodeMerge, nodeLeft, nodeRight, ea, eb);

    //第6步，求RC模型的数值解/求RLC模型的数值解
    if(RLC == false)
        RC_calculation(nodeMerge, nodeLeft, nodeRight, &ea, &eb, L, slope);
    else
        RLC_calculation(nodeMerge, nodeLeft, nodeRight, &ea, &eb, L, slope);

    //printf("/****Test Merge Delay****/\n""Delay before merge: Left: %f, Right: %f, after: Left: %f,  Right: %f\n", nodeLeft->delay,nodeRight->delay, t_a + nodeLeft->delay, t_b + nodeRight->delay);
}


//TODO: 使用RLC的计算公式取代目前使用的naive计算公式
void DME_3D(TreeNode* node, bool RLC){
    if(node){
        //如果已经是根节点，跳过
        if(!node->left_child && !node->right_child){
            return;
        }
        //自下而上递归
        DME_3D(node->left_child, RLC);
        DME_3D(node->right_child, RLC);

        printf("Now update:\n"
               "First: %d, Second: %d, Merge: %d\n",
               node->left_child->ID, node->right_child->ID, node->ID);
        //根据RC模型计算merging point的坐标, 并且更新电容，更新delay
        if(node->left_child->x <= node->right_child->x)
            calculation(node->left_child, node->right_child, node, RLC);
        else
            calculation(node->right_child, node->left_child, node, RLC);

        printf("merge after: x=%.0f, y=%.0f, C: %f\n\n",
               node->x, node->y, node->C);
    }
//
}




#endif //ZST_DME_DME_H
