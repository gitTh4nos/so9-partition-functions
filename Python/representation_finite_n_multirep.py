#!/usr/bin/env python3
"""
Standalone exact finite-n multiplicity for one selected B4/SO(9)
representation.

Usage:

    python3 representation_finite_n.py LEVEL L M R S [options]

All four Dynkin labels are mandatory.  At fixed level n, the program divides
the factorized ZQ directly inside the finite theta-weight support and extracts
only a_[L,M,R,S](n).  It does not need multiplicities of higher
representations.

This is a true single-file program.  It embeds the sparse theta engine,
runtime-compiled pthread C kernel, B4/Weyl utilities, exact bounds, prime/CRT
utilities, checkpointing, Python/Mathematica/JSON output, and optional
verification against the full triangular solver.  It does not import any
other local Python file.
"""

from __future__ import annotations

import argparse
import array
import base64
import ctypes
import concurrent.futures
import hashlib
import itertools
import json
import math
import os
import pickle
import platform
import shutil
import subprocess
import tempfile
import time
import zlib
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, NamedTuple, Optional, Sequence, Tuple


Weight = Tuple[int, int, int, int]
Rep = Tuple[int, int, int, int]

MOD = (1 << 61) - 1
RHO2: Weight = (7, 5, 3, 1)
ETA_STAR: Weight = (4, 0, 0, 0)
INV2 = (MOD + 1) // 2


# ---------------------------------------------------------------------------
# B4 labels, ansatz, and Weyl data
# ---------------------------------------------------------------------------

def ansatz_reps(n: int) -> List[Rep]:
    """Return the same BasisNC B4 ansatz as the supplied dense solver."""
    reps: List[Rep] = []
    for a in range(n + 1):
        for b in range(n // 3 + 1):
            for c in range((n + 2) // 6 + 1):
                for d in range(n // 2 + 1):
                    if (
                        a + 3 * b + 2 * d <= n
                        and a + 3 * b + 3 * d <= n + 1
                        and a + 3 * b + 6 * c + 3 * d <= n + 2
                        and a + 3 * b + 6 * c + 4 * d <= n + 3
                        and a + 3 * b + 6 * c + 5 * d <= n + 7
                    ):
                        reps.append((a, b, c, d))
    return sorted(reps)


def add_weight(a: Weight, b: Weight) -> Weight:
    return tuple(a[i] + b[i] for i in range(4))  # type: ignore[return-value]


def sub_weight(a: Weight, b: Weight) -> Weight:
    return tuple(a[i] - b[i] for i in range(4))  # type: ignore[return-value]


def mu_of_rep(rep: Rep) -> Weight:
    ell, m, r, s = rep
    lam2 = (
        2 * ell + 2 * m + 2 * r + s,
        2 * m + 2 * r + s,
        2 * r + s,
        s,
    )
    return add_weight(lam2, RHO2)


def rep_of_dominant_mu(mu: Weight) -> Optional[Rep]:
    lam = sub_weight(mu, RHO2)
    ell2 = lam[0] - lam[1]
    m2 = lam[1] - lam[2]
    r2 = lam[2] - lam[3]
    s = lam[3]
    if min(ell2, m2, r2, s) < 0:
        return None
    if ell2 % 2 or m2 % 2 or r2 % 2:
        return None
    return (ell2 // 2, m2 // 2, r2 // 2, s)


def permutation_sign(perm: Sequence[int]) -> int:
    inversions = sum(
        perm[i] > perm[j]
        for i in range(len(perm))
        for j in range(i + 1, len(perm))
    )
    return -1 if inversions % 2 else 1


def weyl_rho_shifts() -> List[Tuple[Weight, int]]:
    """Return all 384 signed-permutation images of rho2 and their determinants."""
    result: List[Tuple[Weight, int]] = []
    for perm in itertools.permutations(range(4)):
        p_sign = permutation_sign(perm)
        base = tuple(RHO2[perm[i]] for i in range(4))
        for signs in itertools.product((-1, 1), repeat=4):
            weight = tuple(signs[i] * base[i] for i in range(4))
            sign = p_sign * signs[0] * signs[1] * signs[2] * signs[3]
            result.append((weight, sign))  # type: ignore[arg-type]
    assert len(result) == 384
    assert len({weight for weight, _ in result}) == 384
    return result


WEYL_RHO = weyl_rho_shifts()
WEYL_SHIFT_ARRAY = array.array(
    "i", itertools.chain.from_iterable(weight for weight, _ in WEYL_RHO)
)


def signed_dominant_rep(weight: Weight) -> Tuple[Optional[Rep], int]:
    """
    Reflect a doubled B4 weight into the strict dominant chamber.

    Returns (None,0) on a Weyl wall or when the dominant point is not lambda+rho
    for an integral nonnegative B4 Dynkin label.
    """
    absolute = [abs(value) for value in weight]
    if 0 in absolute or len(set(absolute)) < 4:
        return None, 0
    order = tuple(sorted(range(4), key=lambda i: absolute[i], reverse=True))
    dominant: Weight = tuple(absolute[i] for i in order)  # type: ignore[assignment]
    sign = permutation_sign(order)
    for value in weight:
        sign *= 1 if value > 0 else -1
    return rep_of_dominant_mu(dominant), sign


# ---------------------------------------------------------------------------
# Fixed ZQ weight polynomial
# ---------------------------------------------------------------------------

def zq_weights() -> Dict[Weight, int]:
    """Expand the exact 256-term ZQ product into 153 distinct Laurent weights."""
    zero: Weight = (0, 0, 0, 0)
    z1: Weight = (2, 0, 0, 0)
    z2: Weight = (2, 2, 0, 0)
    z3: Weight = (1, 1, 1, -1)
    z4: Weight = (1, 1, 1, 1)
    factors = (
        (zero, z3),
        (z1, z3),
        (z2, z3),
        (z2, add_weight(z1, z3)),
        (zero, z4),
        (z1, z4),
        (z2, z4),
        (z2, add_weight(z1, z4)),
    )
    denominator: Weight = (12, 8, 4, 0)
    weights: Counter[Weight] = Counter()
    for choices in itertools.product((0, 1), repeat=8):
        exponent = zero
        for factor, choice in zip(factors, choices):
            exponent = add_weight(exponent, factor[choice])
        weights[sub_weight(exponent, denominator)] += 1
    assert len(weights) == 153
    assert sum(weights.values()) == 256
    assert max(weights) == ETA_STAR and weights[ETA_STAR] == 1
    return dict(weights)


ZQ_WEIGHTS = zq_weights()
ZQ_TAIL = tuple((eta, coefficient) for eta, coefficient in ZQ_WEIGHTS.items()
                if eta != ETA_STAR)


# ---------------------------------------------------------------------------
# Runtime-compiled parallel C theta-weight kernel
# ---------------------------------------------------------------------------

_C_SOURCE = r"""
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

static inline uint64_t addmod(uint64_t a, uint64_t b, uint64_t M) {
    uint64_t r = a + b;
    if (r >= M) r -= M;
    return r;
}

static inline uint64_t submod(uint64_t a, uint64_t b, uint64_t M) {
    return a >= b ? a - b : a + M - b;
}

static inline uint64_t mulmod(uint64_t a, uint64_t b, uint64_t M) {
    __uint128_t p = (__uint128_t)a * b;
    const uint64_t mask=UINT64_C(2305843009213693951);
    if (M <= mask && mask-M <= 1024) {
        uint64_t c=(mask+1)-M; /* M=2^61-c. */
        __uint128_t r=(p&mask)+(p>>61)*c;
        r=(r&mask)+(r>>61)*c;
        r=(r&mask)+(r>>61)*c;
        uint64_t value=(uint64_t)r;
        while (value>=M) value-=M;
        return value;
    }
    return (uint64_t)(p % M);
}

static inline size_t tidx(int sector, int ei, int q, int W, int L) {
    return ((size_t)sector * W + ei) * L + q;
}

/*
 * Build three one-variable bivariate tables:
 *   sector 0: G3 = F*C
 *   sector 1: G4 = F*D
 *   sector 2: GB = F*B
 * where F=(x-x^-1)/A in the notation of theta_ab_coeffs.
 */
int build_onevar_tables(int K, int E, uint64_t M, uint64_t* out) {
    const int L = 2*K + 3;
    const int W = 2*E + 1;
    uint64_t* F = (uint64_t*)calloc((size_t)L * W, sizeof(uint64_t));
    if (!F) return 1;
    memset(out, 0, (size_t)3 * W * L * sizeof(uint64_t));
    F[(size_t)0 * W + E] = 1;

    /*
     * H=A/(x-x^-1)
     *  =1+sum_{n>=1}(-1)^n sum_{e=-2n step 2}^{2n} x^e Q^{n(n+1)}
     * and F=1/H, hence F_q=-sum_{j>=1}H_j F_{q-j}.
     */
    for (int q = 1; q < L; q++) {
        uint64_t* dst = F + (size_t)q * W;
        for (int n = 1; n*(n+1) <= q; n++) {
            int j = n*(n+1);
            const uint64_t* src = F + (size_t)(q-j) * W;
            for (int eh = -2*n; eh <= 2*n; eh += 2) {
                int lo = eh < 0 ? -eh : 0;
                int hi = eh > 0 ? W-eh : W;
                for (int ei = lo; ei < hi; ei++) {
                    uint64_t value = src[ei];
                    if (!value) continue;
                    int oi = ei + eh;
                    /* F_q -= (-1)^n * x^eh * F_{q-j}. */
                    if (n & 1) dst[oi] = addmod(dst[oi], value, M);
                    else       dst[oi] = submod(dst[oi], value, M);
                }
            }
        }
    }

    /* Sparse multiplication of F by C, D, and B. */
    for (int qf = 0; qf < L; qf++) {
        const uint64_t* src = F + (size_t)qf * W;
        for (int ei = 0; ei < W; ei++) {
            uint64_t value = src[ei];
            if (!value) continue;

            /* C_0=D_0=1. */
            size_t i3 = tidx(0, ei, qf, W, L);
            size_t i4 = tidx(1, ei, qf, W, L);
            out[i3] = addmod(out[i3], value, M);
            out[i4] = addmod(out[i4], value, M);

            for (int n = 1; qf + n*n < L; n++) {
                int q = qf + n*n;
                int exps[2] = {-2*n, 2*n};
                for (int z = 0; z < 2; z++) {
                    int oi = ei + exps[z];
                    if (oi < 0 || oi >= W) continue;
                    i3 = tidx(0, oi, q, W, L);
                    i4 = tidx(1, oi, q, W, L);
                    out[i3] = addmod(out[i3], value, M);
                    if (n & 1) out[i4] = submod(out[i4], value, M);
                    else       out[i4] = addmod(out[i4], value, M);
                }
            }

            for (int n = 0; qf + n*(n+1) < L; n++) {
                int q = qf + n*(n+1);
                int e = 2*n + 1;
                int exps[2] = {-e, e};
                for (int z = 0; z < 2; z++) {
                    int oi = ei + exps[z];
                    if (oi < 0 || oi >= W) continue;
                    size_t ib = tidx(2, oi, q, W, L);
                    out[ib] = addmod(out[ib], value, M);
                }
            }
        }
    }
    free(F);
    return 0;
}

static uint64_t plan_hash4(const int32_t* weight) {
    uint64_t h=UINT64_C(0x9e3779b97f4a7c15);
    for (int i=0;i<4;i++) {
        uint64_t x=(uint64_t)(uint32_t)weight[i]+UINT64_C(0x9e3779b9);
        x ^= x >> 16; x *= UINT64_C(0x85ebca6b);
        x ^= x >> 13; x *= UINT64_C(0xc2b2ae35);
        x ^= x >> 16;
        h ^= x+UINT64_C(0x9e3779b97f4a7c15)+(h<<6)+(h>>2);
    }
    return h;
}

/*
 * Build the 384 Weyl-shifted theta requests for each representation and
 * deduplicate them in compiled code. out_ids has N*S entries; out_unique has
 * capacity 4*N*S. The returned unique count selects its valid prefix.
 */
int required_theta_plan_batch(
    const int32_t* reps, int N, const int32_t* shifts, int S,
    int32_t* out_unique, int32_t* out_ids, int32_t* out_unique_count
) {
    size_t max_items=(size_t)N*(size_t)S;
    size_t hash_size=1;
    while (hash_size < 2*max_items) hash_size <<= 1;
    int32_t* hash_values=(int32_t*)calloc(hash_size,sizeof(int32_t));
    int32_t* hash_weights=(int32_t*)malloc(4*hash_size*sizeof(int32_t));
    if (!hash_values || !hash_weights) {
        free(hash_values); free(hash_weights); return 1;
    }
    size_t mask=hash_size-1;
    int32_t unique_count=0;
    for (int i=0;i<N;i++) {
        const int32_t* p=reps+4*(size_t)i;
        int32_t gamma[4]={
            2*p[0]+2*p[1]+2*p[2]+p[3]+11,
            2*p[1]+2*p[2]+p[3]+5,
            2*p[2]+p[3]+3,
            p[3]+1
        };
        for (int j=0;j<S;j++) {
            const int32_t* shift=shifts+4*(size_t)j;
            int32_t weight[4]={
                gamma[0]-shift[0],gamma[1]-shift[1],
                gamma[2]-shift[2],gamma[3]-shift[3]
            };
            size_t slot=(size_t)plan_hash4(weight)&mask;
            for (;;) {
                int32_t stored=hash_values[slot];
                if (!stored) {
                    int32_t index=unique_count++;
                    memcpy(out_unique+4*(size_t)index,weight,4*sizeof(int32_t));
                    memcpy(hash_weights+4*slot,weight,4*sizeof(int32_t));
                    hash_values[slot]=index+1;
                    out_ids[(size_t)i*S+j]=index;
                    break;
                }
                const int32_t* known=hash_weights+4*slot;
                if (known[0]==weight[0] && known[1]==weight[1] &&
                    known[2]==weight[2] && known[3]==weight[3]) {
                    out_ids[(size_t)i*S+j]=stored-1;
                    break;
                }
                slot=(slot+1)&mask;
            }
        }
    }
    *out_unique_count=unique_count;
    free(hash_values); free(hash_weights);
    return 0;
}

static uint64_t plan_hash2(int32_t a, int32_t b) {
    int32_t weight[4]={a,b,0,0};
    return plan_hash4(weight);
}

/* Deduplicate the two canonical exponent pairs of every unique weight. */
int required_pair_plan_batch(
    const int32_t* weights, int U,
    int32_t* out_pairs, int32_t* out_pair_ids, int32_t* out_pair_count
) {
    size_t max_pairs=2*(size_t)U;
    size_t hash_size=1;
    while (hash_size < 2*max_pairs) hash_size <<= 1;
    int32_t* hash_values=(int32_t*)calloc(hash_size,sizeof(int32_t));
    int32_t* hash_pairs=(int32_t*)malloc(2*hash_size*sizeof(int32_t));
    if (!hash_values || !hash_pairs) {
        free(hash_values); free(hash_pairs); return 1;
    }
    size_t mask=hash_size-1;
    int32_t pair_count=0;
    for (int i=0;i<U;i++) {
        const int32_t* weight=weights+4*(size_t)i;
        for (int side=0;side<2;side++) {
            int32_t a=weight[2*side],b=weight[2*side+1];
            if (a>b) { int32_t tmp=a; a=b; b=tmp; }
            size_t slot=(size_t)plan_hash2(a,b)&mask;
            for (;;) {
                int32_t stored=hash_values[slot];
                if (!stored) {
                    int32_t index=pair_count++;
                    out_pairs[2*(size_t)index]=a;
                    out_pairs[2*(size_t)index+1]=b;
                    hash_pairs[2*slot]=a; hash_pairs[2*slot+1]=b;
                    hash_values[slot]=index+1;
                    out_pair_ids[2*(size_t)i+side]=index;
                    break;
                }
                if (hash_pairs[2*slot]==a && hash_pairs[2*slot+1]==b) {
                    out_pair_ids[2*(size_t)i+side]=stored-1;
                    break;
                }
                slot=(slot+1)&mask;
            }
        }
    }
    *out_pair_count=pair_count;
    free(hash_values); free(hash_pairs);
    return 0;
}

typedef struct {
    const uint64_t* table;
    const int32_t* pairs;
    uint64_t* out;
    int L, W, E, U, start, end;
    uint64_t M;
} PairArgs;

static void* pair_worker(void* ptr) {
    PairArgs* a = (PairArgs*)ptr;
    for (int u = a->start; u < a->end; u++) {
        int e1 = a->pairs[2*u] + a->E;
        int e2 = a->pairs[2*u+1] + a->E;
        if (e1 < 0 || e1 >= a->W || e2 < 0 || e2 >= a->W) continue;
        for (int sector = 0; sector < 3; sector++) {
            const uint64_t* left =
                a->table + ((size_t)sector*a->W + e1)*a->L;
            const uint64_t* right =
                a->table + ((size_t)sector*a->W + e2)*a->L;
            uint64_t* dst =
                a->out + ((size_t)sector*a->U + u)*a->L;
            for (int q = 0; q < a->L; q++) {
                uint64_t acc = 0;
                for (int i = 0; i <= q; i++) {
                    if (!left[i] || !right[q-i]) continue;
                    acc = addmod(acc, mulmod(left[i], right[q-i], a->M), a->M);
                }
                dst[q] = acc;
            }
        }
    }
    return NULL;
}

int pair_series_batch(
    const uint64_t* table, int L, int W, int E,
    const int32_t* pairs, int U, uint64_t M,
    uint64_t* out, int n_threads
) {
    memset(out, 0, (size_t)3*U*L*sizeof(uint64_t));
    if (U == 0) return 0;
    if (n_threads < 1) n_threads = 1;
    if (n_threads > U) n_threads = U;
    pthread_t* threads = (pthread_t*)malloc((size_t)n_threads*sizeof(pthread_t));
    PairArgs* args = (PairArgs*)malloc((size_t)n_threads*sizeof(PairArgs));
    if (!threads || !args) { free(threads); free(args); return 1; }
    int chunk = (U + n_threads - 1)/n_threads;
    for (int t = 0; t < n_threads; t++) {
        int start = t*chunk, end = start+chunk;
        if (end > U) end = U;
        args[t] = (PairArgs){table,pairs,out,L,W,E,U,start,end,M};
        pthread_create(&threads[t],NULL,pair_worker,&args[t]);
    }
    for (int t = 0; t < n_threads; t++) pthread_join(threads[t],NULL);
    free(threads); free(args);
    return 0;
}

typedef struct {
    const uint64_t* pair_series;
    const int32_t* pair_ids;
    uint64_t* out;
    int K, L, U, N, start, end;
    uint64_t M;
} ThetaArgs;

static uint64_t pair_product_coeff(
    const uint64_t* pair_series, int sector, int U, int L,
    int p1, int p2, int target, uint64_t M
) {
    const uint64_t* a = pair_series + ((size_t)sector*U+p1)*L;
    const uint64_t* b = pair_series + ((size_t)sector*U+p2)*L;
    uint64_t acc = 0;
    for (int i = 0; i <= target; i++) {
        if (!a[i] || !b[target-i]) continue;
        acc = addmod(acc,mulmod(a[i],b[target-i],M),M);
    }
    return acc;
}

static void* theta_worker(void* ptr) {
    ThetaArgs* a = (ThetaArgs*)ptr;
    uint64_t inv2 = (a->M + 1)/2;
    int ta = 2*a->K + 1;
    int tb = 2*a->K;
    for (int i = a->start; i < a->end; i++) {
        int p1 = a->pair_ids[2*i], p2 = a->pair_ids[2*i+1];
        uint64_t g3 = pair_product_coeff(a->pair_series,0,a->U,a->L,p1,p2,ta,a->M);
        uint64_t g4 = pair_product_coeff(a->pair_series,1,a->U,a->L,p1,p2,ta,a->M);
        uint64_t gb = pair_product_coeff(a->pair_series,2,a->U,a->L,p1,p2,tb,a->M);
        uint64_t value = addmod(submod(g3,g4,a->M),gb,a->M);
        a->out[i] = mulmod(value,inv2,a->M);
    }
    return NULL;
}

int theta_from_pairs_batch(
    const uint64_t* pair_series, int K, int L, int U,
    const int32_t* pair_ids, int N, uint64_t M,
    uint64_t* out, int n_threads
) {
    if (N == 0) return 0;
    if (n_threads < 1) n_threads = 1;
    if (n_threads > N) n_threads = N;
    pthread_t* threads = (pthread_t*)malloc((size_t)n_threads*sizeof(pthread_t));
    ThetaArgs* args = (ThetaArgs*)malloc((size_t)n_threads*sizeof(ThetaArgs));
    if (!threads || !args) { free(threads); free(args); return 1; }
    int chunk = (N+n_threads-1)/n_threads;
    for (int t = 0; t < n_threads; t++) {
        int start=t*chunk,end=start+chunk;
        if (end>N) end=N;
        args[t]=(ThetaArgs){pair_series,pair_ids,out,K,L,U,N,start,end,M};
        pthread_create(&threads[t],NULL,theta_worker,&args[t]);
    }
    for (int t=0;t<n_threads;t++) pthread_join(threads[t],NULL);
    free(threads); free(args);
    return 0;
}

typedef struct {
    const uint64_t* theta_values;
    const int32_t* weyl_ids;
    const int32_t* weyl_signs;
    uint64_t* out;
    int N, start, end;
    uint64_t M;
} WeylArgs;

static void* weyl_worker(void* ptr) {
    WeylArgs* a = (WeylArgs*)ptr;
    for (int i = a->start; i < a->end; i++) {
        uint64_t acc = 0;
        const int32_t* ids = a->weyl_ids + (size_t)i*384;
        for (int j = 0; j < 384; j++) {
            uint64_t value = a->theta_values[ids[j]];
            if (a->weyl_signs[j] > 0) acc = addmod(acc,value,a->M);
            else                      acc = submod(acc,value,a->M);
        }
        a->out[i] = acc;
    }
    return NULL;
}

int weyl_sum_batch(
    const uint64_t* theta_values, const int32_t* weyl_ids,
    const int32_t* weyl_signs, int N, uint64_t M,
    uint64_t* out, int n_threads
) {
    if (N == 0) return 0;
    if (n_threads < 1) n_threads = 1;
    if (n_threads > N) n_threads = N;
    pthread_t* threads = (pthread_t*)malloc((size_t)n_threads*sizeof(pthread_t));
    WeylArgs* args = (WeylArgs*)malloc((size_t)n_threads*sizeof(WeylArgs));
    if (!threads || !args) { free(threads); free(args); return 1; }
    int chunk = (N+n_threads-1)/n_threads;
    for (int t=0;t<n_threads;t++) {
        int start=t*chunk,end=start+chunk;
        if (end>N) end=N;
        args[t]=(WeylArgs){theta_values,weyl_ids,weyl_signs,out,N,start,end,M};
        pthread_create(&threads[t],NULL,weyl_worker,&args[t]);
    }
    for (int t=0;t<n_threads;t++) pthread_join(threads[t],NULL);
    free(threads); free(args);
    return 0;
}

typedef struct {
    int N;
    int base;
    uint64_t M;
    int32_t* reps;
    uint64_t* recovered;
    uint64_t* hash_keys;
    int32_t* hash_values;
    size_t hash_size;
} RecurrenceContext;

static uint64_t mix64(uint64_t x) {
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

static uint64_t rep_key(int l, int m, int r, int s, int base) {
    uint64_t key = (uint64_t)l;
    key = key*(uint64_t)base + (uint64_t)m;
    key = key*(uint64_t)base + (uint64_t)r;
    key = key*(uint64_t)base + (uint64_t)s;
    return key + 1;
}

static int rep_lookup(const RecurrenceContext* ctx, int l, int m, int r, int s) {
    if (l < 0 || m < 0 || r < 0 || s < 0 ||
        l >= ctx->base || m >= ctx->base || r >= ctx->base || s >= ctx->base)
        return -1;
    uint64_t key = rep_key(l,m,r,s,ctx->base);
    size_t mask = ctx->hash_size-1;
    size_t slot = (size_t)mix64(key) & mask;
    while (ctx->hash_keys[slot]) {
        if (ctx->hash_keys[slot] == key) return ctx->hash_values[slot];
        slot = (slot+1)&mask;
    }
    return -1;
}

void* recurrence_context_create(
    const int32_t* reps, int N, int base, uint64_t M
) {
    RecurrenceContext* ctx = (RecurrenceContext*)calloc(1,sizeof(RecurrenceContext));
    if (!ctx) return NULL;
    ctx->N=N; ctx->base=base; ctx->M=M;
    ctx->reps=(int32_t*)malloc((size_t)4*N*sizeof(int32_t));
    ctx->recovered=(uint64_t*)calloc((size_t)N,sizeof(uint64_t));
    ctx->hash_size=1;
    while (ctx->hash_size < (size_t)2*N) ctx->hash_size <<= 1;
    ctx->hash_keys=(uint64_t*)calloc(ctx->hash_size,sizeof(uint64_t));
    ctx->hash_values=(int32_t*)malloc(ctx->hash_size*sizeof(int32_t));
    if (!ctx->reps || !ctx->recovered || !ctx->hash_keys || !ctx->hash_values) {
        free(ctx->reps); free(ctx->recovered); free(ctx->hash_keys);
        free(ctx->hash_values); free(ctx); return NULL;
    }
    memcpy(ctx->reps,reps,(size_t)4*N*sizeof(int32_t));
    size_t mask=ctx->hash_size-1;
    for (int i=0;i<N;i++) {
        const int32_t* p=ctx->reps+4*(size_t)i;
        uint64_t key=rep_key(p[0],p[1],p[2],p[3],base);
        size_t slot=(size_t)mix64(key)&mask;
        while (ctx->hash_keys[slot]) slot=(slot+1)&mask;
        ctx->hash_keys[slot]=key;
        ctx->hash_values[slot]=i;
    }
    return ctx;
}

static int reflected_rep(
    const int alpha_in[4], int* l, int* m, int* r, int* s, int* sign
) {
    int a[4];
    int sg=1;
    for (int i=0;i<4;i++) {
        if (alpha_in[i]==0) return 0;
        if (alpha_in[i]<0) { a[i]=-alpha_in[i]; sg=-sg; }
        else a[i]=alpha_in[i];
    }
    for (int i=0;i<4;i++)
        for (int j=i+1;j<4;j++) {
            if (a[i]==a[j]) return 0;
            if (a[i]<a[j]) {
                int tmp=a[i]; a[i]=a[j]; a[j]=tmp; sg=-sg;
            }
        }
    int lam0=a[0]-7,lam1=a[1]-5,lam2=a[2]-3,lam3=a[3]-1;
    int dl=lam0-lam1,dm=lam1-lam2,dr=lam2-lam3;
    if (dl<0 || dm<0 || dr<0 || lam3<0 ||
        (dl&1) || (dm&1) || (dr&1)) return 0;
    *l=dl/2; *m=dm/2; *r=dr/2; *s=lam3; *sign=sg;
    return 1;
}

int recurrence_recover_batch(
    void* opaque, int start, int count, const uint64_t* rhs,
    const int32_t* zq_weights, const int32_t* zq_coeffs, int Z
) {
    RecurrenceContext* ctx=(RecurrenceContext*)opaque;
    if (!ctx || start<0 || count<0 || start+count>ctx->N) return 1;
    for (int local=0;local<count;local++) {
        int index=start+local;
        const int32_t* p=ctx->reps+4*(size_t)index;
        int l0=p[0],m0=p[1],r0=p[2],s0=p[3];
        int mu[4]={
            2*l0+2*m0+2*r0+s0+7,
            2*m0+2*r0+s0+5,
            2*r0+s0+3,
            s0+1
        };
        int gamma[4]={mu[0]+4,mu[1],mu[2],mu[3]};
        uint64_t correction=0;
        for (int z=0;z<Z;z++) {
            const int32_t* eta=zq_weights+4*(size_t)z;
            int alpha[4]={
                gamma[0]-eta[0],gamma[1]-eta[1],
                gamma[2]-eta[2],gamma[3]-eta[3]
            };
            int l,m,r,s,sg;
            if (!reflected_rep(alpha,&l,&m,&r,&s,&sg)) continue;
            int neighbor=rep_lookup(ctx,l,m,r,s);
            if (neighbor<0) continue;
            if (neighbor>=index) return 2;
            uint64_t term=mulmod(
                (uint64_t)zq_coeffs[z],ctx->recovered[neighbor],ctx->M
            );
            if (sg>0) correction=addmod(correction,term,ctx->M);
            else      correction=submod(correction,term,ctx->M);
        }
        ctx->recovered[index]=submod(rhs[local],correction,ctx->M);
    }
    return 0;
}

int recurrence_copy_values(void* opaque, uint64_t* out) {
    RecurrenceContext* ctx=(RecurrenceContext*)opaque;
    if (!ctx) return 1;
    memcpy(out,ctx->recovered,(size_t)ctx->N*sizeof(uint64_t));
    return 0;
}

void recurrence_context_destroy(void* opaque) {
    RecurrenceContext* ctx=(RecurrenceContext*)opaque;
    if (!ctx) return;
    free(ctx->reps); free(ctx->recovered); free(ctx->hash_keys);
    free(ctx->hash_values); free(ctx);
}
"""


def _compile_c_kernel(cache_dir: Optional[Path] = None) -> ctypes.CDLL:
    compiler = shutil.which("clang") or shutil.which("cc") or shutil.which("gcc")
    if compiler is None:
        raise RuntimeError("A C compiler (clang, cc, or gcc) is required")

    digest = hashlib.sha256(_C_SOURCE.encode("utf-8")).hexdigest()[:16]
    suffix = ".dylib" if platform.system() == "Darwin" else ".so"
    if cache_dir is None:
        cache_dir = Path(tempfile.gettempdir()) / "pf_weyl_sparse"
    cache_dir.mkdir(parents=True, exist_ok=True)
    source = cache_dir / f"weyl_sparse_{digest}.c"
    library = cache_dir / f"weyl_sparse_{digest}{suffix}"
    if not library.exists():
        source.write_text(_C_SOURCE, encoding="utf-8")
        command = [
            compiler,
            "-O3",
            "-std=c11",
            "-shared",
            "-fPIC",
            "-pthread",
            str(source),
            "-o",
            str(library),
        ]
        subprocess.run(command, check=True)
    lib = ctypes.CDLL(str(library))

    u64p = ctypes.POINTER(ctypes.c_uint64)
    i32p = ctypes.POINTER(ctypes.c_int32)
    lib.build_onevar_tables.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_uint64, u64p
    ]
    lib.build_onevar_tables.restype = ctypes.c_int
    lib.required_theta_plan_batch.argtypes = [
        i32p, ctypes.c_int, i32p, ctypes.c_int, i32p, i32p, i32p,
    ]
    lib.required_theta_plan_batch.restype = ctypes.c_int
    lib.required_pair_plan_batch.argtypes = [
        i32p, ctypes.c_int, i32p, i32p, i32p,
    ]
    lib.required_pair_plan_batch.restype = ctypes.c_int
    lib.pair_series_batch.argtypes = [
        u64p, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        i32p, ctypes.c_int, ctypes.c_uint64, u64p, ctypes.c_int,
    ]
    lib.pair_series_batch.restype = ctypes.c_int
    lib.theta_from_pairs_batch.argtypes = [
        u64p, ctypes.c_int, ctypes.c_int, ctypes.c_int,
        i32p, ctypes.c_int, ctypes.c_uint64, u64p, ctypes.c_int,
    ]
    lib.theta_from_pairs_batch.restype = ctypes.c_int
    lib.weyl_sum_batch.argtypes = [
        u64p, i32p, i32p, ctypes.c_int, ctypes.c_uint64, u64p, ctypes.c_int,
    ]
    lib.weyl_sum_batch.restype = ctypes.c_int
    lib.recurrence_context_create.argtypes = [
        i32p, ctypes.c_int, ctypes.c_int, ctypes.c_uint64,
    ]
    lib.recurrence_context_create.restype = ctypes.c_void_p
    lib.recurrence_recover_batch.argtypes = [
        ctypes.c_void_p, ctypes.c_int, ctypes.c_int, u64p, i32p, i32p,
        ctypes.c_int,
    ]
    lib.recurrence_recover_batch.restype = ctypes.c_int
    lib.recurrence_copy_values.argtypes = [ctypes.c_void_p, u64p]
    lib.recurrence_copy_values.restype = ctypes.c_int
    lib.recurrence_context_destroy.argtypes = [ctypes.c_void_p]
    lib.recurrence_context_destroy.restype = None
    return lib


def _u64_pointer(values: array.array) -> ctypes.POINTER(ctypes.c_uint64):
    return ctypes.cast(
        (ctypes.c_uint64 * len(values)).from_buffer(values),
        ctypes.POINTER(ctypes.c_uint64),
    )


def _i32_pointer(values: array.array) -> ctypes.POINTER(ctypes.c_int32):
    return ctypes.cast(
        (ctypes.c_int32 * len(values)).from_buffer(values),
        ctypes.POINTER(ctypes.c_int32),
    )


class SparseThetaEngine:
    """Parallel sparse theta-weight evaluator for one fixed q order."""

    def __init__(
        self,
        order: int,
        exponent_bound: int,
        *,
        modulus: int = MOD,
        threads: int = 1,
        c_cache_dir: Optional[Path] = None,
    ):
        self.order = order
        self.L = 2 * order + 3
        self.E = exponent_bound
        self.W = 2 * exponent_bound + 1
        self.modulus = modulus
        self.threads = max(1, threads)
        self.lib = _compile_c_kernel(c_cache_dir)

        table_size = 3 * self.W * self.L
        self.table = array.array("Q", [0]) * table_size
        status = self.lib.build_onevar_tables(
            order, self.E, modulus, _u64_pointer(self.table)
        )
        if status:
            raise MemoryError("build_onevar_tables failed")

    @staticmethod
    def _canonical_pair(a: int, b: int) -> Tuple[int, int]:
        return (a, b) if a <= b else (b, a)

    def evaluate(self, weights: Sequence[Weight]) -> Dict[Weight, int]:
        """Evaluate unique theta coefficients for one in-memory weight batch."""
        if not weights:
            return {}
        unique_weights = list(dict.fromkeys(weights))
        flat_weights = array.array(
            "i", itertools.chain.from_iterable(unique_weights)
        )
        output = self.evaluate_flat(flat_weights, len(unique_weights))
        return dict(zip(unique_weights, map(int, output)))

    def evaluate_flat(
        self,
        flat_weights: array.array,
        unique_count: int,
    ) -> array.array:
        """Evaluate a compiled flat array of unique four-component weights."""
        if unique_count == 0:
            return array.array("Q")
        if len(flat_weights) < 4 * unique_count:
            raise ValueError("flat theta-weight array has the wrong length")
        if any(abs(flat_weights[index]) > self.E for index in range(4 * unique_count)):
            raise ValueError(f"theta weight exceeds exponent bound E={self.E}")

        pair_to_id: Dict[Tuple[int, int], int] = {}
        pair_ids = array.array("i")
        for index in range(unique_count):
            offset = 4 * index
            left = self._canonical_pair(
                flat_weights[offset], flat_weights[offset + 1]
            )
            right = self._canonical_pair(
                flat_weights[offset + 2], flat_weights[offset + 3]
            )
            for pair in (left, right):
                if pair not in pair_to_id:
                    pair_to_id[pair] = len(pair_to_id)
            pair_ids.append(pair_to_id[left])
            pair_ids.append(pair_to_id[right])

        pairs_by_id: List[Tuple[int, int]] = [(-1, -1)] * len(pair_to_id)
        for pair, index in pair_to_id.items():
            pairs_by_id[index] = pair
        pairs = array.array("i", itertools.chain.from_iterable(pairs_by_id))
        return self.evaluate_pair_plan(
            pairs, pair_ids, unique_count, len(pairs_by_id)
        )

    def evaluate_pair_plan(
        self,
        pairs: array.array,
        pair_ids: array.array,
        unique_count: int,
        pair_count: int,
    ) -> array.array:
        """Evaluate theta weights from a pre-deduplicated exponent-pair plan."""
        pair_series = array.array("Q", [0]) * (3 * pair_count * self.L)
        status = self.lib.pair_series_batch(
            _u64_pointer(self.table),
            self.L,
            self.W,
            self.E,
            _i32_pointer(pairs),
            pair_count,
            self.modulus,
            _u64_pointer(pair_series),
            self.threads,
        )
        if status:
            raise MemoryError("pair_series_batch failed")

        output = array.array("Q", [0]) * unique_count
        status = self.lib.theta_from_pairs_batch(
            _u64_pointer(pair_series),
            self.order,
            self.L,
            pair_count,
            _i32_pointer(pair_ids),
            unique_count,
            self.modulus,
            _u64_pointer(output),
            self.threads,
        )
        if status:
            raise RuntimeError("theta_from_pairs_batch failed")
        return output

    def weyl_sums(
        self,
        theta_values: Sequence[int],
        weyl_ids: array.array,
        n_reps: int,
    ) -> List[int]:
        values = array.array("Q", theta_values)
        signs = array.array("i", (sign for _, sign in WEYL_RHO))
        output = array.array("Q", [0]) * n_reps
        status = self.lib.weyl_sum_batch(
            _u64_pointer(values),
            _i32_pointer(weyl_ids),
            _i32_pointer(signs),
            n_reps,
            self.modulus,
            _u64_pointer(output),
            self.threads,
        )
        if status:
            raise RuntimeError("weyl_sum_batch failed")
        return list(map(int, output))


# ---------------------------------------------------------------------------
# Triangular multiplicity solver
# ---------------------------------------------------------------------------

class TriangularRecurrenceContext:
    """Compiled ordered solve of ZQ * multiplicities = Weyl numerator."""

    def __init__(
        self,
        lib: ctypes.CDLL,
        descending_reps: Sequence[Rep],
        order: int,
        modulus: int,
    ):
        self.lib = lib
        self.count = len(descending_reps)
        self._context: Optional[int] = None
        self.reps = array.array(
            "i", itertools.chain.from_iterable(descending_reps)
        )
        self.zq_weights = array.array(
            "i", itertools.chain.from_iterable(eta for eta, _ in ZQ_TAIL)
        )
        self.zq_coeffs = array.array(
            "i", (coefficient for _, coefficient in ZQ_TAIL)
        )
        context = self.lib.recurrence_context_create(
            _i32_pointer(self.reps), self.count, order + 1, modulus
        )
        if not context:
            raise MemoryError("recurrence_context_create failed")
        self._context = context

    def recover(self, start: int, rhs_values: Sequence[int]) -> None:
        rhs = array.array("Q", rhs_values)
        status = self.lib.recurrence_recover_batch(
            self._context,
            start,
            len(rhs),
            _u64_pointer(rhs),
            _i32_pointer(self.zq_weights),
            _i32_pointer(self.zq_coeffs),
            len(ZQ_TAIL),
        )
        if status == 2:
            raise RuntimeError(
                "Triangular order violation in compiled recurrence"
            )
        if status:
            raise RuntimeError("recurrence_recover_batch failed")

    def values(self) -> List[int]:
        output = array.array("Q", [0]) * self.count
        status = self.lib.recurrence_copy_values(
            self._context, _u64_pointer(output)
        )
        if status:
            raise RuntimeError("recurrence_copy_values failed")
        return list(map(int, output))

    def close(self) -> None:
        if self._context is not None:
            self.lib.recurrence_context_destroy(self._context)
            self._context = None

    def __enter__(self) -> "TriangularRecurrenceContext":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def chunks(values: Sequence[Rep], size: int) -> Iterator[Sequence[Rep]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def required_theta_plan(reps: Sequence[Rep]) -> Tuple[List[Weight], array.array]:
    """Reference Python implementation retained for tests."""
    requested: List[Weight] = []
    request_index: Dict[Weight, int] = {}
    weyl_ids = array.array("i")
    for rep in reps:
        gamma = add_weight(mu_of_rep(rep), ETA_STAR)
        for shift, _ in WEYL_RHO:
            weight = sub_weight(gamma, shift)
            index = request_index.get(weight)
            if index is None:
                index = len(requested)
                request_index[weight] = index
                requested.append(weight)
            weyl_ids.append(index)
    return requested, weyl_ids


class CompiledThetaPlan(NamedTuple):
    pairs: array.array
    pair_ids: array.array
    weyl_ids: array.array
    unique_count: int
    pair_count: int


def required_theta_plan_compiled(
    lib: ctypes.CDLL,
    reps: Sequence[Rep],
) -> CompiledThetaPlan:
    """Compiled generation and deduplication of all Weyl-shifted requests."""
    n_reps = len(reps)
    max_items = n_reps * len(WEYL_RHO)
    flat_reps = array.array("i", itertools.chain.from_iterable(reps))
    unique_weights = array.array("i", [0]) * (4 * max_items)
    weyl_ids = array.array("i", [0]) * max_items
    unique_count = ctypes.c_int32()
    status = lib.required_theta_plan_batch(
        _i32_pointer(flat_reps),
        n_reps,
        _i32_pointer(WEYL_SHIFT_ARRAY),
        len(WEYL_RHO),
        _i32_pointer(unique_weights),
        _i32_pointer(weyl_ids),
        ctypes.byref(unique_count),
    )
    if status:
        raise MemoryError("required_theta_plan_batch failed")
    count = int(unique_count.value)

    pairs = array.array("i", [0]) * (4 * count)
    pair_ids = array.array("i", [0]) * (2 * count)
    pair_count = ctypes.c_int32()
    status = lib.required_pair_plan_batch(
        _i32_pointer(unique_weights),
        count,
        _i32_pointer(pairs),
        _i32_pointer(pair_ids),
        ctypes.byref(pair_count),
    )
    if status:
        raise MemoryError("required_pair_plan_batch failed")
    return CompiledThetaPlan(
        pairs=pairs,
        pair_ids=pair_ids,
        weyl_ids=weyl_ids,
        unique_count=count,
        pair_count=int(pair_count.value),
    )


def exponent_bound_for_order(order: int) -> int:
    # Highest ansatz mu1 is at most 2*order+7; gamma and the Weyl shift add 11.
    return 2 * order + 20


def solve_order_moduli(
    order: int,
    moduli: Sequence[int],
    *,
    threads: int = 1,
    rep_batch_size: int = 1000,
    progress: bool = True,
    c_cache_dir: Optional[Path] = None,
) -> List[Dict[Rep, int]]:
    """
    Compute one q order for several moduli while sharing every batch plan.

    rep_batch_size controls peak pair-series memory.  Larger batches reuse more
    theta weights/pairs but require more RAM.
    """
    if not moduli:
        raise ValueError("at least one modulus is required")
    all_reps = ansatz_reps(order)
    descending = sorted(all_reps, key=mu_of_rep, reverse=True)
    engines = [
        SparseThetaEngine(
            order + 1,
            exponent_bound_for_order(order + 1),
            modulus=modulus,
            threads=threads,
            c_cache_dir=c_cache_dir,
        )
        for modulus in moduli
    ]
    recurrences = [
        TriangularRecurrenceContext(
            engine.lib, descending, order, modulus
        )
        for engine, modulus in zip(engines, moduli)
    ]

    started = time.perf_counter()
    done = 0
    try:
        for rep_batch in chunks(descending, rep_batch_size):
            plan = required_theta_plan_compiled(
                engines[0].lib, rep_batch
            )
            for engine, recurrence in zip(engines, recurrences):
                theta_values = engine.evaluate_pair_plan(
                    plan.pairs,
                    plan.pair_ids,
                    plan.unique_count,
                    plan.pair_count,
                )
                rhs_values = engine.weyl_sums(
                    theta_values, plan.weyl_ids, len(rep_batch)
                )
                recurrence.recover(done, rhs_values)
            done += len(rep_batch)

            if progress:
                elapsed = time.perf_counter() - started
                print(
                    f"q={order}: {done:,}/{len(descending):,} reps "
                    f"({done/elapsed:,.1f} reps/s), "
                    f"theta_weights={plan.unique_count:,}, "
                    f"moduli={len(moduli)}",
                    flush=True,
                )
        recovered_by_modulus = [
            recurrence.values() for recurrence in recurrences
        ]
    finally:
        for recurrence in recurrences:
            recurrence.close()

    return [
        {
            rep: value
            for rep, value in zip(descending, recovered_values)
            if value
        }
        for recovered_values in recovered_by_modulus
    ]


def solve_order(
    order: int,
    *,
    threads: int = 1,
    rep_batch_size: int = 1000,
    modulus: int = MOD,
    progress: bool = True,
    c_cache_dir: Optional[Path] = None,
) -> Dict[Rep, int]:
    """Compute all nonzero B4 multiplicities at one q order modulo modulus."""
    return solve_order_moduli(
        order,
        (modulus,),
        threads=threads,
        rep_batch_size=rep_batch_size,
        progress=progress,
        c_cache_dir=c_cache_dir,
    )[0]


# ---------------------------------------------------------------------------
# Exact bounds, multiple primes, and CRT
# ---------------------------------------------------------------------------

def _integer_convolution(a: Sequence[int], b: Sequence[int], length: int) -> List[int]:
    out = [0] * length
    for i, left in enumerate(a):
        if not left:
            continue
        for j in range(min(len(b), length - i)):
            right = b[j]
            if right:
                out[i + j] += left * right
    return out


@lru_cache(maxsize=None)
def coefficient_upper_bound(order: int) -> int:
    """
    Rigorous upper bound from the unrefined coefficient Z_k(1).

    Z_k(1)=sum_lambda a_{lambda,k} dim(lambda), so positivity and
    dim(lambda)>=1 imply a_{lambda,k} <= Z_k(1).
    """
    theta_order = order + 1
    length = 2 * theta_order + 3

    # H(1)=1+sum (-1)^n (2n+1) Q^(n(n+1)); F=1/H.
    h = [0] * length
    h[0] = 1
    n = 1
    while n * (n + 1) < length:
        h[n * (n + 1)] += (-1 if n & 1 else 1) * (2 * n + 1)
        n += 1
    f = [0] * length
    f[0] = 1
    for q in range(1, length):
        f[q] = -sum(h[j] * f[q - j] for j in range(1, q + 1))

    c = [0] * length
    d = [0] * length
    b = [0] * length
    c[0] = d[0] = 1
    n = 1
    while n * n < length:
        c[n * n] += 2
        d[n * n] += -2 if n & 1 else 2
        n += 1
    n = 0
    while n * (n + 1) < length:
        b[n * (n + 1)] += 2
        n += 1

    g3 = _integer_convolution(f, c, length)
    g4 = _integer_convolution(f, d, length)
    gb = _integer_convolution(f, b, length)

    def fourth_power(series: Sequence[int]) -> List[int]:
        square = _integer_convolution(series, series, length)
        return _integer_convolution(square, square, length)

    g3_four = fourth_power(g3)
    g4_four = fourth_power(g4)
    gb_four = fourth_power(gb)
    theta_at_one = (
        g3_four[2 * theta_order + 1]
        - g4_four[2 * theta_order + 1]
        + gb_four[2 * theta_order]
    ) // 2
    if theta_at_one % 256:
        raise ArithmeticError("unrefined theta coefficient is not divisible by ZQ(1)=256")
    bound = theta_at_one // 256
    if bound < 0:
        raise ArithmeticError("negative unrefined state bound")
    return bound


def is_prime_64(value: int) -> bool:
    """Deterministic Miller-Rabin test for unsigned 64-bit integers."""
    if value < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small_primes:
        if value % prime == 0:
            return value == prime
    odd_part = value - 1
    twos = 0
    while odd_part % 2 == 0:
        twos += 1
        odd_part //= 2
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(twos - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


@lru_cache(maxsize=None)
def reconstruction_primes(count: int) -> Tuple[int, ...]:
    """Return distinct primes descending from the 61-bit Mersenne prime."""
    if count < 1:
        raise ValueError("at least one reconstruction prime is required")
    primes: List[int] = []
    candidate = MOD
    while len(primes) < count:
        if is_prime_64(candidate):
            primes.append(candidate)
        candidate -= 2
    return tuple(primes)


def primes_for_bound(bound: int) -> Tuple[int, ...]:
    product = 1
    count = 0
    while product <= bound:
        count += 1
        product *= reconstruction_primes(count)[-1]
    return reconstruction_primes(max(1, count))


_B4_DIMENSION_DENOMINATOR = math.prod(RHO2)
for _dimension_i in range(4):
    for _dimension_j in range(_dimension_i + 1, 4):
        _B4_DIMENSION_DENOMINATOR *= (
            RHO2[_dimension_i] ** 2 - RHO2[_dimension_j] ** 2
        )


def b4_dimension(rep: Rep) -> int:
    """Exact Weyl dimension of a B4 representation."""
    mu = mu_of_rep(rep)
    numerator = math.prod(mu)
    for i in range(4):
        for j in range(i + 1, 4):
            numerator *= mu[i] * mu[i] - mu[j] * mu[j]
    quotient, remainder = divmod(numerator, _B4_DIMENSION_DENOMINATOR)
    if remainder:
        raise ArithmeticError(f"nonintegral B4 dimension for {rep}")
    return quotient


def unrefined_state_sum(coefficients: Mapping[Rep, int]) -> int:
    return sum(
        coefficient * b4_dimension(rep)
        for rep, coefficient in coefficients.items()
    )


def crt_reconstruct(
    residue_maps: Sequence[Mapping[Rep, int]],
    primes: Sequence[int],
) -> Tuple[Dict[Rep, int], int]:
    """Incremental CRT for sparse, nonnegative coefficient dictionaries."""
    if not residue_maps or len(residue_maps) != len(primes):
        raise ValueError("one residue dictionary is required for each prime")
    values: Dict[Rep, int] = dict(residue_maps[0])
    product = primes[0]
    for residue_map, prime in zip(residue_maps[1:], primes[1:]):
        inverse = pow(product % prime, -1, prime)
        existing = list(values)
        for rep in existing:
            residue = residue_map.get(rep, 0)
            correction = ((residue - values[rep] % prime) * inverse) % prime
            values[rep] += product * correction
        for rep, residue in residue_map.items():
            if rep not in values:
                values[rep] = product * ((residue * inverse) % prime)
        product *= prime
    return {rep: value for rep, value in values.items() if value}, product


def solve_order_exact(
    order: int,
    *,
    threads: int = 1,
    rep_batch_size: int = 1000,
    progress: bool = True,
    verify_prime: bool = False,
    prime_workers: int = 1,
    c_cache_dir: Optional[Path] = None,
) -> Tuple[Dict[Rep, int], int, Tuple[int, ...]]:
    """Compute exact nonnegative multiplicities using an automatic CRT basis."""
    bound = coefficient_upper_bound(order)
    primes = primes_for_bound(bound)
    if progress:
        print(
            f"q={order}: rigorous coefficient bound={bound} "
            f"({bound.bit_length()} bits), "
            f"at most {len(primes)} CRT primes",
            flush=True,
        )
    residues: List[Dict[Rep, int]] = []
    exact: Dict[Rep, int] = {}
    used_primes: Tuple[int, ...] = ()
    prime_workers = max(1, min(prime_workers, len(primes)))
    if prime_workers > 1:
        # Ensure the shared library exists before worker threads try to load it.
        _compile_c_kernel(c_cache_dir)

    solved = False
    next_prime = 0
    while next_prime < len(primes) and not solved:
        group = primes[next_prime : next_prime + prime_workers]
        worker_threads = max(1, threads // len(group))
        if progress:
            print(
                f"q={order}: evaluating primes {next_prime+1}-"
                f"{next_prime+len(group)}/{len(primes)} "
                f"with {len(group)} prime worker(s)",
                flush=True,
            )

        def solve_prime(prime: int) -> Dict[Rep, int]:
            return solve_order(
                order,
                threads=worker_threads,
                rep_batch_size=rep_batch_size,
                modulus=prime,
                progress=progress and len(group) == 1,
                c_cache_dir=c_cache_dir,
            )

        if len(group) == 1:
            group_residues = [solve_prime(group[0])]
        else:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=len(group)
            ) as executor:
                group_residues = list(executor.map(solve_prime, group))

        for residue in group_residues:
            residues.append(residue)
            prime_index = len(residues)
            used_primes = primes[:prime_index]
            exact, product = crt_reconstruct(residues, used_primes)
            state_sum = unrefined_state_sum(exact)
            if state_sum > bound:
                raise ArithmeticError(
                    "CRT representative exceeds the unrefined state count"
                )
            if progress:
                print(
                    f"q={order}: dimension checksum after {prime_index} "
                    f"prime(s): {state_sum}/{bound}",
                    flush=True,
                )
            if state_sum == bound:
                solved = True
                break
        next_prime += len(group)

    if not solved:
        raise ArithmeticError(
            "CRT product exhausted the rigorous bound but dimension checksum failed"
        )

    if exact and max(exact.values()) > bound:
        raise ArithmeticError("reconstructed coefficient exceeds the unrefined bound")

    if verify_prime:
        check_prime = reconstruction_primes(len(used_primes) + 1)[-1]
        check = solve_order(
            order,
            threads=threads,
            rep_batch_size=rep_batch_size,
            modulus=check_prime,
            progress=progress,
            c_cache_dir=c_cache_dir,
        )
        for rep in set(exact) | set(check):
            if exact.get(rep, 0) % check_prime != check.get(rep, 0):
                raise ArithmeticError(f"CRT verification failed for {rep}")
        if progress:
            print(f"q={order}: independent prime verification passed", flush=True)
    return exact, bound, used_primes


# ---------------------------------------------------------------------------
# Validation and I/O
# ---------------------------------------------------------------------------

def validate_against(
    result: Mapping[Rep, int],
    expected: Mapping[Rep, int],
    *,
    modulus: Optional[int] = MOD,
) -> List[Tuple[Rep, int, int]]:
    reps = set(result) | set(expected)
    mismatches = []
    for rep in sorted(reps):
        got = result.get(rep, 0)
        want = expected.get(rep, 0)
        if modulus is not None:
            got %= modulus
            want %= modulus
        if got != want:
            mismatches.append((rep, got, want))
    return mismatches


def save_pickle(path: Path, results: Mapping[int, Mapping[Rep, int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as stream:
        serializable = {
            order: dict(coefficients) for order, coefficients in results.items()
        }
        pickle.dump(serializable, stream, protocol=pickle.HIGHEST_PROTOCOL)


def _mx_escape_string(value: str) -> bytes:
    """Port of the original solver's Mathematica String encoder."""
    out = bytearray()
    for character in value:
        code = ord(character)
        if code < 256:
            out.append(code)
        else:
            utf16 = character.encode("utf-16-be")
            for index in range(0, len(utf16), 2):
                unit = int.from_bytes(utf16[index : index + 2], "big")
                out.extend(f"\\:{unit:04X}".encode("ascii"))
    return bytes(out)


def mathematica_compress(value: str, level: int = 9) -> str:
    """Return a Mathematica Compress-compatible payload for a String."""
    escaped = _mx_escape_string(value)
    encoded_string = b"!boR" + b"S" + len(escaped).to_bytes(4, "little") + escaped
    return "1:" + base64.b64encode(zlib.compress(encoded_string, level)).decode("ascii")


def mathematica_association(results: Mapping[int, Mapping[Rep, int]]) -> str:
    """Serialize as <|n -> <|{l,m,r,s} -> multiplicity, ...|>, ...|>."""
    levels = []
    for order in sorted(results):
        entries = []
        for rep, coefficient in sorted(results[order].items()):
            labels = ",".join(map(str, rep))
            entries.append(f"{{{labels}}}->{coefficient}")
        levels.append(f"{order}-><|" + ",".join(entries) + "|>")
    return "<|" + ",".join(levels) + "|>"


def save_mathematica_association(
    path: Path,
    results: Mapping[int, Mapping[Rep, int]],
    *,
    compressed: bool = True,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    association = mathematica_association(results)
    if compressed:
        payload = mathematica_compress(association)
        content = f'ToExpression[Uncompress["{payload}"]]\n'
    else:
        content = association + "\n"
    path.write_text(content, encoding="ascii")


def cache_path(cache_dir: Path, order: int, exact: bool) -> Path:
    kind = "exact" if exact else "modM"
    return cache_dir / f"level_{order:04d}_{kind}.pkl"


def load_cached_level(
    cache_dir: Path, order: int, exact: bool
) -> Optional[Dict[Rep, int]]:
    path = cache_path(cache_dir, order, exact)
    if not path.exists():
        return None
    with path.open("rb") as stream:
        cached = pickle.load(stream)
    return dict(cached[order])


def save_cached_level(
    cache_dir: Path,
    order: int,
    coefficients: Mapping[Rep, int],
    exact: bool,
) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    save_pickle(cache_path(cache_dir, order, exact), {order: coefficients})


# ---------------------------------------------------------------------------
# Direct finite division for one selected representation
# ---------------------------------------------------------------------------

DELTAS: Tuple[Weight, ...] = tuple(
    (1, e2, e3, e4)
    for e2 in (-1, 1)
    for e3 in (-1, 1)
    for e4 in (-1, 1)
)


def _combine_rules(rules: Iterable[Tuple[Weight, int]]) -> Dict[Weight, int]:
    result: Counter[Weight] = Counter()
    for weight, coefficient in rules:
        result[weight] += coefficient
    return {
        weight: coefficient
        for weight, coefficient in result.items()
        if coefficient
    }


def _divide_by_one_plus_monomial(
    old: Mapping[Weight, int],
    delta: Weight,
    first_coordinate_max: int,
) -> Dict[Weight, int]:
    """Truncated exact division by (1+y^(-delta))."""
    layers: Dict[int, Dict[Tuple[int, int, int], int]] = {}
    for (first, second, third, fourth), coefficient in old.items():
        layers.setdefault(first, {})[(second, third, fourth)] = coefficient

    if not layers:
        return {}

    d2, d3, d4 = delta[1:]
    previous: Dict[Tuple[int, int, int], int] = {}
    output: Dict[Weight, int] = {}
    first_min = min(layers)

    for first in range(first_min, first_coordinate_max + 1):
        current: Counter[Tuple[int, int, int]] = Counter(
            layers.get(first, {})
        )
        for (second, third, fourth), coefficient in previous.items():
            shifted = (second + d2, third + d3, fourth + d4)
            current[shifted] -= coefficient
        previous = {
            triple: coefficient
            for triple, coefficient in current.items()
            if coefficient
        }
        for triple, coefficient in previous.items():
            output[(first, *triple)] = coefficient

    return output


def representation_projector_weights(
    level: int,
    representation: Rep,
) -> Dict[Weight, int]:
    """Return p_nu such that a_rep(n)=Sum_nu p_nu Theta_n[nu]."""
    if level < 0:
        raise ValueError("level must be nonnegative")
    if any(label < 0 for label in representation):
        raise ValueError("B4 Dynkin labels must be nonnegative")

    theta_support = 2 * level + 4
    gamma = add_weight(mu_of_rep(representation), ETA_STAR)

    projector = _combine_rules(
        (
            tuple(gamma[i] - image[i] for i in range(4)),
            sign,
        )
        for image, sign in WEYL_RHO
    )
    projector = {
        weight: coefficient
        for weight, coefficient in projector.items()
        if weight[0] <= theta_support
    }

    for delta in DELTAS:
        projector = _divide_by_one_plus_monomial(
            projector, delta, theta_support
        )

    projector = {
        weight: coefficient
        for weight, coefficient in projector.items()
        if all(abs(component) <= theta_support for component in weight)
    }

    return _combine_rules(
        (
            tuple(sorted(map(abs, weight), reverse=True)),
            coefficient,
        )
        for weight, coefficient in projector.items()
    )


def representation_mod_prime(
    level: int,
    projector: Mapping[Weight, int],
    prime: int,
    *,
    threads: int,
) -> int:
    """Evaluate one selected-representation projector modulo one prime."""
    engine = SparseThetaEngine(
        level + 1,
        exponent_bound_for_order(level + 1),
        modulus=prime,
        threads=threads,
    )
    theta_values = engine.evaluate(list(projector))
    result = 0
    for weight, coefficient in projector.items():
        result = (result + coefficient * theta_values[weight]) % prime
    return result


def crt_pair(
    value: int,
    product: int,
    residue: int,
    prime: int,
) -> Tuple[int, int]:
    """Extend x=value (mod product) with x=residue (mod prime)."""
    correction = (
        (residue - value % prime) * pow(product % prime, -1, prime)
    ) % prime
    return value + product * correction, product * prime


def representation_exact(
    level: int,
    representation: Rep,
    *,
    threads: int = 1,
    prime_workers: int = 1,
    progress: bool = True,
    projector: Mapping[Weight, int] | None = None,
) -> Tuple[int, int, Tuple[int, ...], int]:
    """Compute and rigorously reconstruct a_[representation](level)."""
    started = time.perf_counter()
    if projector is None:
        projector = representation_projector_weights(level, representation)
    build_elapsed = time.perf_counter() - started

    if progress:
        max_coefficient = max(map(abs, projector.values()), default=0)
        print(
            f"n={level}, rep={representation}: projector has "
            f"{len(projector):,} theta weights; "
            f"max |projector coefficient|={max_coefficient}; "
            f"built in {build_elapsed:.3f}s",
            flush=True,
        )

    bound = coefficient_upper_bound(level)
    primes = primes_for_bound(bound)
    prime_workers = max(1, min(prime_workers, len(primes)))
    worker_threads = max(1, threads // prime_workers)
    if prime_workers > 1:
        _compile_c_kernel(None)

    def evaluate_prime(prime: int) -> int:
        return representation_mod_prime(
            level,
            projector,
            prime,
            threads=worker_threads,
        )

    if prime_workers == 1:
        residues = [evaluate_prime(prime) for prime in primes]
    else:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=prime_workers
        ) as executor:
            residues = list(executor.map(evaluate_prime, primes))

    value = 0
    product = 1
    for index, (prime, residue) in enumerate(zip(primes, residues), 1):
        value, product = crt_pair(value, product, residue, prime)
        if progress:
            print(
                f"n={level}, rep={representation}: prime "
                f"{index}/{len(primes)} -> residue {residue}; "
                f"CRT product {product.bit_length()} bits",
                flush=True,
            )

    if product <= bound:
        raise ArithmeticError("CRT product does not exceed the rigorous bound")
    if not 0 <= value <= bound:
        raise ArithmeticError(
            f"reconstructed multiplicity {value} lies outside "
            f"0 <= a <= {bound}"
        )
    return value, bound, primes, len(projector)


def representation_cache_path(
    cache_dir: Path,
    level: int,
    representation: Rep,
) -> Path:
    labels = "_".join(map(str, representation))
    return cache_dir / f"rep_{labels}_level_{level:04d}_exact.pkl"


def load_representation_cache(
    cache_dir: Path,
    level: int,
    representation: Rep,
) -> Optional[dict]:
    path = representation_cache_path(cache_dir, level, representation)
    if not path.exists():
        return None
    with path.open("rb") as stream:
        record = pickle.load(stream)
    if record.get("level") != level:
        raise ValueError(f"invalid level in representation cache {path}")
    cached_representation = tuple(record.get("representation", ()))
    if cached_representation != representation:
        raise ValueError(
            f"cache representation {cached_representation} does not match "
            f"requested representation {representation} in {path}"
        )
    return record


def save_representation_cache(cache_dir: Path, record: Mapping) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    representation = tuple(record["representation"])
    if len(representation) != 4:
        raise ValueError("cache representation must contain four Dynkin labels")
    path = representation_cache_path(
        cache_dir,
        int(record["level"]),
        representation,  # type: ignore[arg-type]
    )
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as stream:
        pickle.dump(dict(record), stream, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)
# ---------------------------------------------------------------------------
# Log-file helper
# ---------------------------------------------------------------------------

import logging  # add to the imports at the top of the file if not already there

def _setup_logger(log_path: Optional[Path], quiet: bool) -> logging.Logger:
    """
    Return a logger that always writes to *log_path* (if given) and also
    writes to stdout unless *quiet* is True.

    Each line is timestamped:
        2024-03-15 12:34:56  INFO  n=5, rep=(1,0,0,0): completed in 3.217s
    """
    logger = logging.getLogger("rep_finite_n")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()          # avoid duplicate handlers on re-entry

    formatter = logging.Formatter(
        fmt="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    if not quiet:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    # If no handlers were added (quiet=True and no log_path) use NullHandler
    if not logger.handlers:
        logger.addHandler(logging.NullHandler())

    return logger
# ---------------------------------------------------------------------------
# Multi-representation file parser  (insert just above main())
# ---------------------------------------------------------------------------

def parse_representations_file(path: Path) -> List[Rep]:
    """
    Parse a list of B4 representations from a text file.

    Each non-blank, non-comment line must contain exactly four non-negative
    integers (Dynkin labels l, m, r, s) separated by whitespace or commas.

    Example file::

        # l  m  r  s
        0  0  0  0
        1  0  0  0
        0, 1, 0, 0
    """
    representations: List[Rep] = []
    with path.open("r", encoding="utf-8") as stream:
        for line_num, raw in enumerate(stream, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) != 4:
                raise ValueError(
                    f"{path}:{line_num}: expected 4 integers (l m r s), "
                    f"got {len(parts)}: {line!r}"
                )
            try:
                rep: Rep = tuple(int(p) for p in parts)  # type: ignore[assignment]
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_num}: cannot parse integers from {line!r}"
                ) from exc
            if any(label < 0 for label in rep):
                raise ValueError(
                    f"{path}:{line_num}: Dynkin labels must be "
                    f"nonneg: {rep}"
                )
            representations.append(rep)
    if not representations:
        raise ValueError(f"no valid representations found in {path}")
    return representations


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Standalone exact multiplicity of one or more B4 representations "
            "from the finite theta formula.  Supply l m r s positionally for a "
            "single representation, or pass --rep-file for multiple."
        )
    )
    parser.add_argument("level", type=int, help="stored q-order n")
    parser.add_argument(
        "labels",
        nargs="*",
        type=int,
        metavar="LABEL",
        help="four B4 Dynkin labels l m r s (required unless --rep-file is given)",
    )
    parser.add_argument(
        "--rep-file",
        type=Path,
        metavar="FILE",
        default=None,
        help=(
            "text file with one representation per line "
            "(four integers l m r s, space- or comma-separated); "
            "when given, any positional LABEL arguments are ignored"
        ),
    )
    parser.add_argument(
        "--per-rep-dir",
        type=Path,
        metavar="DIR",
        default=None,
        help="write a separate JSON result file per representation into this directory",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="threads used by the sparse theta backend (default: 1)",
    )
    parser.add_argument(
        "--prime-workers",
        type=int,
        default=1,
        help="CRT primes evaluated concurrently (default: 1)",
    )
    parser.add_argument(
        "--up-to",
        action="store_true",
        help="compute every level from --start-level through level",
    )
    parser.add_argument(
        "--start-level",
        type=int,
        default=1,
        help="first level used with --up-to (default: 1)",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=None,
        help="resumable per-level cache; filenames include all four labels",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="write {n: {(l,m,r,s): multiplicity}} as a Python pickle",
    )
    parser.add_argument(
        "--mathematica-output",
        type=Path,
        default=None,
        help="write the nested Mathematica Association used by the full solver",
    )
    parser.add_argument(
        "--mathematica-plain",
        action="store_true",
        help="do not Compress the Mathematica Association",
    )
    parser.add_argument(
        "--verify-full",
        action="store_true",
        help="run the embedded full solver and compare the selected result",
    )
    parser.add_argument(
        "--json-output",
        type=Path,
        default=None,
        help="write the result and certification metadata as JSON",
    )
    # ── new ──────────────────────────────────────────────────────────────────
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        metavar="FILE",
        help=(
            "append timestamped completion records to this file; "
            "each line records the representation, multiplicity, and "
            "wall-clock time taken"
        ),
    )
    # ─────────────────────────────────────────────────────────────────────────
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    # ── logger ────────────────────────────────────────────────────────────────
    logger = _setup_logger(args.log_file, args.quiet)

    # ── build the list of representations ─────────────────────────────────────
    representations: List[Rep]
    if args.rep_file is not None:
        try:
            representations = parse_representations_file(args.rep_file)
        except (OSError, ValueError) as exc:
            parser.error(str(exc))
            return
        if args.labels and not args.quiet:
            print(
                "Warning: positional LABEL arguments are ignored "
                "when --rep-file is given.",
                flush=True,
            )
    elif len(args.labels) == 4:
        if any(v < 0 for v in args.labels):
            parser.error("all four B4 Dynkin labels must be nonnegative")
        representations = [tuple(args.labels)]  # type: ignore[list-item]
    else:
        parser.error(
            f"expected 4 Dynkin labels (l m r s) but got {len(args.labels)}; "
            "use --rep-file to supply multiple representations"
        )
        return

    # ── shared validation ──────────────────────────────────────────────────────
    if args.level < 0:
        parser.error("level must be nonnegative")
    if args.start_level < 0:
        parser.error("--start-level must be nonnegative")
    if args.up_to and args.start_level > args.level:
        parser.error("--start-level cannot exceed level")

    started = time.perf_counter()
    levels: List[int] = (
        list(range(args.start_level, args.level + 1))
        if args.up_to
        else [args.level]
    )

    logger.info(
        "starting: %d representation(s), level(s) %d-%d",
        len(representations),
        levels[0],
        levels[-1],
    )

    aggregate: Dict[int, Dict[Rep, int]] = {level: {} for level in levels}
    all_metadata: List[dict] = []

    # ── outer loop: one pass per representation ────────────────────────────────
    for rep_index, representation in enumerate(representations, 1):
        rep_started = time.perf_counter()                        # ← new

        if not args.quiet and len(representations) > 1:
            print(
                f"\n── representation {rep_index}/{len(representations)}: "
                f"{representation} ──",
                flush=True,
            )

        logger.info(                                             # ← new
            "representation %d/%d %s: starting",
            rep_index,
            len(representations),
            representation,
        )

        projector_buckets: Dict[int, List[Tuple[Weight, int]]] = {}
        active_projector: Optional[Dict[Weight, int]] = None
        active_support = -1

        if args.up_to:
            master_projector = representation_projector_weights(
                args.level, representation
            )
            for weight, coefficient in master_projector.items():
                projector_buckets.setdefault(max(weight), []).append(
                    (weight, coefficient)
                )
            active_projector = {}

        rep_metadata: List[dict] = []

        # ── inner loop: one pass per level ────────────────────────────────────
        for level in levels:
            level_started = time.perf_counter()

            if active_projector is not None:
                support = 2 * level + 4
                for coordinate_max in range(active_support + 1, support + 1):
                    active_projector.update(
                        projector_buckets.get(coordinate_max, ())
                    )
                active_support = support

            cached = (
                load_representation_cache(args.cache_dir, level, representation)
                if args.cache_dir
                else None
            )

            if cached is None:
                multiplicity, bound, primes, weight_count = representation_exact(
                    level,
                    representation,
                    threads=max(1, args.threads),
                    prime_workers=max(1, args.prime_workers),
                    progress=not args.quiet,
                    projector=active_projector,
                )
                record: dict = {
                    "level": level,
                    "representation": list(representation),
                    "multiplicity": multiplicity,
                    "rigorous_upper_bound": bound,
                    "crt_primes": list(primes),
                    "theta_weight_count": weight_count,
                }
                if args.cache_dir:
                    save_representation_cache(args.cache_dir, record)
            else:
                record = cached
                multiplicity = int(record["multiplicity"])
                if not args.quiet:
                    print(
                        f"n={level}, rep={representation}: loaded from cache",
                        flush=True,
                    )

            if args.verify_full:
                full, _, _ = solve_order_exact(
                    level,
                    threads=max(1, args.threads),
                    progress=not args.quiet,
                    prime_workers=max(1, args.prime_workers),
                )
                expected = full.get(representation, 0)
                if multiplicity != expected:
                    raise ArithmeticError(
                        f"selected result {multiplicity} != full solver result "
                        f"{expected} for rep={representation}, n={level}"
                    )
                if not args.quiet:
                    print(
                        f"n={level}, rep={representation}: "
                        "full-solver verification passed",
                        flush=True,
                    )

            if multiplicity:
                aggregate[level][representation] = multiplicity

            level_elapsed = time.perf_counter() - level_started   # ← new
            record = dict(record)
            record["elapsed_seconds"] = level_elapsed
            rep_metadata.append(record)

            if not args.quiet:
                print(
                    f"n={level}: a_{representation}={multiplicity} "
                    f"({level_elapsed:.3f}s)",
                    flush=True,
                )
        # ── end inner level loop ───────────────────────────────────────────────

        rep_elapsed = time.perf_counter() - rep_started           # ← new
        all_metadata.extend(rep_metadata)

        # ── log completion of this representation ─────────────────────────────
        # Build a compact multiplicity summary for the log line, e.g.
        #   n=5->3  or  n=1->0, n=2->1, n=3->2, n=4->3, n=5->6
        mult_summary = ", ".join(
            f"n={rec['level']}->{rec['multiplicity']}"
            for rec in rep_metadata
        )
        logger.info(                                              # ← new
            "representation %d/%d %s: COMPLETED in %.3fs | %s",
            rep_index,
            len(representations),
            representation,
            rep_elapsed,
            mult_summary,
        )

        # Optional per-representation JSON output
        if args.per_rep_dir is not None:
            args.per_rep_dir.mkdir(parents=True, exist_ok=True)
            rep_label = "_".join(map(str, representation))
            per_path = args.per_rep_dir / f"rep_{rep_label}_levels.json"
            per_path.write_text(
                json.dumps(
                    {
                        "representation": list(representation),
                        "results": rep_metadata,
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            if not args.quiet:
                print(f"saved {per_path}", flush=True)
    # ── end outer representation loop ─────────────────────────────────────────

    elapsed = time.perf_counter() - started

    # ── aggregate file outputs ─────────────────────────────────────────────────
    if args.output:
        save_pickle(args.output, aggregate)
        if not args.quiet:
            print(f"saved Python aggregate {args.output}", flush=True)

    if args.mathematica_output:
        save_mathematica_association(
            args.mathematica_output,
            aggregate,
            compressed=not args.mathematica_plain,
        )
        if not args.quiet:
            print(
                f"saved Mathematica Association {args.mathematica_output}",
                flush=True,
            )

    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(
            json.dumps(
                {
                    "representations": [list(rep) for rep in representations],
                    "levels": levels,
                    "results": all_metadata,
                    "elapsed_seconds": elapsed,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
        if not args.quiet:
            print(f"saved JSON output {args.json_output}", flush=True)

    # ── final summary ─────────────────────────────────────────────────────────
    logger.info(                                                  # ← new
        "all done: %d representation(s), levels %d-%d, total %.3fs",
        len(representations),
        levels[0],
        levels[-1],
        elapsed,
    )

    if len(representations) > 1 or args.up_to:
        print(
            f"completed {len(representations)} representation(s), "
            f"levels {levels[0]}-{levels[-1]} in {elapsed:.3f}s"
        )
    else:
        print(all_metadata[0]["multiplicity"])


if __name__ == "__main__":
    main()