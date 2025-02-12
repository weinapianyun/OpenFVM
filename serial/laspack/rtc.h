/****************************************************************************/
/*                                  rtc.h                                   */
/****************************************************************************/
/*                                                                          */
/* Residual Termination Control                                             */
/*                                                                          */
/* Copyright (C) 1992-1996 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef RTC_H
#define RTC_H

#include "lastypes.h"
#include "vector.h"
#include "itersolv.h"
#include "copyrght.h"

/* identifiers for iteration methods */
// 迭代方法的标识符
typedef enum {
    /* classical iterative methods */
	GSIterId,
    JacobiIterId,
    SORForwIterId,
    SORBackwIterId,
    SSORIterId,

    /* semi-iterative methods */
    ChebyshevIterId,

    /* CG and CG-like methods */
    CGIterId,
    CGNIterId,
    GMRESIterId,
    BiCGIterId,
    QMRIterId,
    CGSIterId,
    BiCGSTABIterId,

    /* multigrid and multigrid based methods */
    // 多重网格系列算法
    MGIterId,
    NestedMGIterId,
    MGPCGIterId,
    BPXPCGIterId
} IterIdType;

typedef void (*RTCAuxProcType)(int, double, double, IterIdType);

void SetRTCAccuracy(double Eps); // 设置求解终止残差(RTS)
void SetRTCAuxProc(RTCAuxProcType AuxProc);
Boolean RTCResult(int Iter, double rNorm, double bNorm, IterIdType IterId);
int GetLastNoIter(void); // 获取最后一次迭代后的总迭代次数
double GetLastAccuracy(void); // 获取在迭代最后达到的精度

#endif /* RTC_H */
