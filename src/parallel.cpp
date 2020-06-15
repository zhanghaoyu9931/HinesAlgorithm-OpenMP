#include "omp.h"
#include <cstdio>
#include <cassert>
#include <cstring>
#include <algorithm>

// #define DEBUG

#define HASH_TABLE_SIZE (1 << 6)
#define HASH(i) (i & (HASH_TABLE_SIZE - 1))
omp_lock_t hash_lock[HASH_TABLE_SIZE];
double startTime, endTime;
double t1, t2, t3, t4;

int N;
double *u, *l, *d, *rhs;
int *p;

int *n_chil;      // n_chil[i] is #children of node i.
int *n_chil_done; // how many child nodes are done in backward sweep stage.
bool *isCrit;     // if node i is critical, i.e. shared by multiple threads?
int *next;        // next[i] is the next sibling of node i, 0 means no more sibling.
int *size;        // size[i] is size of the subtree with root node i,
                  // valid only for branching nodes and the first node in each branch.
int *tidAssigned; // tidAssigned[i] is the tid assigned to node i,
                  // valid only for leaf nodes.

struct
{
    int *size; // point to the same array as the global 'size'.
    bool operator()(int i1, int i2)
    {
        return size[i1] < size[i2];
    }
} myCompare;

// Assign thread @tid to leaf nodes between [@snode,@enode)
void AssignLeafNodes(int tid, int snode, int enode)
{
#ifdef DEBUG
    printf("AssignLeafNodes(%d,%d,%d)\n", tid, snode, enode);
#endif // DEBUG
    for (int i = snode; i < enode; i++)
        if (n_chil[i] == 0)
            tidAssigned[i] = tid;
}

// Assign nodes of a subtree to threads:
// Make sure that size[@root] is set beforehand.
// params:
// @root: root node of this subtree, make sure it has multiple children.
// @cutOff: if any branch is smaller than @cutOff in size, then assign
//          it to a single thread.
// @stid: starting thread id to be assigned.
// return:
// @etid: ending thread id, i.e. last assigned tid plus 1.
// side effect:
// set size[i], tidAssigned[i], isCrit[i] for each node i in this subtree.
void AssignNodes(int root, int cutOff, int stid, int *etid)
{
#ifdef DEBUG
    printf("AssignNodes(%d,%d,%d,&)\n", root, cutOff, stid);
    assert(n_chil[root] > 1);
#endif // DEBUG

    int i, nexti, k, r;
    int grouping, tid;
    int *child = new int[n_chil[root]];

    // get child; get size[i] for each child node i of @root
    for (nexti = next[i = root + 1], k = 0; nexti > 0; nexti = next[i = nexti])
    {
        child[k++] = i;
        size[i] = nexti - i;
    }
    child[k++] = i;
    size[i] = root + size[root] - i;

    // sort child in increasing order in branch size
    std::sort(child, child + n_chil[root], myCompare);

    // assign nodes as even as possible
    grouping = 0;
    tid = stid;
    for (k = 0; k < n_chil[root]; k++) // for each branch after @root
    {
        i = child[k];
        if (size[i] <= cutOff) // this branch is too small, let's assign it to
                               // a single thread, possibly grouping together
                               // with other small branches.
        {
            AssignLeafNodes(tid, i, i + size[i]);
            grouping += size[i];
            if (grouping > cutOff) // use another thread
            {
                grouping = 0;
                tid++;
            }
        }
        else // this branch is large, consideration of multiple threads should be made.
        {
            for (r = i; n_chil[r] == 1; r++) // walk along non-branching nodes
                ;
            if (n_chil[r] == 0) // we reach a leaf node
                tidAssigned[r] = tid++;
            else // we find next branching point, i.e. a new subtree with root node r
            {
                size[r] = size[i] - (r - i);
                if (size[r] <= cutOff) // the new subtree is too small,
                                       // let's assign it to a single thread.
                    AssignLeafNodes(tid++, r, r + size[r]);
                else // the new subtree is large, we need recursive assignment.
                    AssignNodes(r, cutOff, tid, &tid);
            }
        }
    }
    *etid = tid;
    if (tid > stid + 1) // more than one thread will update node @root
        isCrit[root] = true;

    delete[] child;
}

// Recursive forward sweep
// params:
// @root: root node of this subtree, make sure it has multiple children.
void ForwardUpdate(int root)
{
#ifdef DEBUG
    assert(n_chil[root] > 1);
#endif // DEBUG

    int i, k;
    for (k = root + 1; k != 0; k = next[k])
    {
#pragma omp task private(i) firstprivate(k)
        {
#ifdef DEBUG
            printf("new task. root=%d. k=%d. tid=%d.\n", root, k, omp_get_thread_num());
#endif // DEBUG
            for (i = k; n_chil[i] == 1; i++)
            {
                rhs[i] -= l[i] * rhs[p[i]];
                rhs[i] /= d[i];
            }
            rhs[i] -= l[i] * rhs[p[i]];
            rhs[i] /= d[i];
            if (n_chil[i] > 1)
                ForwardUpdate(i);
        }
    }
}

void HinesAlgo()
{
    bool flag;
    int i, j, k, tid, n_procs, cutOff, nt;
    double factor, tmp1, tmp2;

    n_chil = new int[N];
    n_chil_done = new int[N];
    isCrit = new bool[N];
    next = new int[N];
    myCompare.size = size = new int[N];
    tidAssigned = new int[N];
    memset(n_chil, 0, N * sizeof(int));
    memset(n_chil_done, 0, N * sizeof(int));
    memset(isCrit, 0, N * sizeof(bool));
    memset(next, 0, N * sizeof(int));
    memset(size, 0, N * sizeof(int));
    memset(tidAssigned, 0, N * sizeof(int));

    for (i = 0; i < HASH_TABLE_SIZE; i++)
        omp_init_lock(&hash_lock[i]);

    n_procs = omp_get_num_procs();
    cutOff = N / (n_procs * 2);
    printf("num of threads available: %d\n", n_procs);
    printf("cutOff set as: %d\n", cutOff);

    t1 = omp_get_wtime();

    // -----------------------------------------------------------------
    // Setup n_chil, next, size, isCrit, tidAssigned.

    for (i = N - 1; i > 0; i--)
    {
        j = p[i];
        n_chil[j]++;
        if (j != i - 1) // new branch found
        {
            next[i] = next[j + 1];
            next[j + 1] = i;
        }
    }

    for (i = 0; n_chil[i] == 1; i++) // find the first branching point
        ;
    if (n_chil[i] == 0) // the whole tree is non-branching, we should apply serial algo
    {
        tidAssigned[i] = 0;
        nt = 1;
    }
    else // node i is the first branching point
    {
        size[i] = N - i;
        AssignNodes(i, cutOff, 0, &nt);
    }
    printf("num of threads used in backward sweep: %d\n", nt);

#ifdef DEBUG
    for (i = 0; i < N; i++)
        printf("size[%d]=%d, tidAssigned[%d]=%d, isCrit[%d]=%d\n",
               i, size[i], i, tidAssigned[i], i, isCrit[i]);
#endif // DEBUG

    t2 = omp_get_wtime();

    // -----------------------------------------------------------------
    // Backward sweep

#pragma omp parallel private(flag, i, j, k, tid, factor, tmp1, tmp2) num_threads(nt)
    {
        tid = omp_get_thread_num();
        for (k = 0; k < N; k++)
        {
            if (n_chil[k] != 0 || tidAssigned[k] != tid)
                continue;
            i = k;
            j = p[i];
            while (j >= 0) // "parent" of the root is -1
            {
                // update d[j], rhs[j] and n_chil_done[j]
                factor = u[i] / d[i];
                tmp1 = factor * l[i];
                tmp2 = factor * rhs[i];

                if (isCrit[j])
                    omp_set_lock(&hash_lock[HASH(j)]);
                d[j] -= tmp1;
                rhs[j] -= tmp2;
                n_chil_done[j]++;
                flag = (n_chil[j] == n_chil_done[j]);
                if (isCrit[j])
                    omp_unset_lock(&hash_lock[HASH(j)]);

                if (!flag) // this branch finishes
                    break;
                else // update i,j
                    j = p[i = j];
            }
        }
    }

    t3 = omp_get_wtime();

    // -----------------------------------------------------------------
    // Forward sweep

    rhs[0] /= d[0];
#pragma omp parallel num_threads(n_procs)
    {
#pragma omp single
        {
            for (i = 0; n_chil[i] == 1; i++)
            {
                if (i != 0)
                {
                    rhs[i] -= l[i] * rhs[p[i]];
                    rhs[i] /= d[i];
                }
            }
            if (i != 0)
            {
                rhs[i] -= l[i] * rhs[p[i]];
                rhs[i] /= d[i];
            }
            if (n_chil[i] > 1)
                ForwardUpdate(i);
        }
    }

    t4 = omp_get_wtime();

    delete[] n_chil;
    delete[] n_chil_done;
    delete[] isCrit;
    delete[] next;
    delete[] size;
    delete[] tidAssigned;
}

int main(int argc, char const *argv[])
{
    int i, idx;
    FILE *fin, *fout;

    // read data from file
    assert(argc == 3);
    fin = fopen(argv[1], "r");
    fscanf(fin, "%d", &N);
    u = new double[N];
    l = new double[N];
    d = new double[N];
    rhs = new double[N];
    p = new int[N];
    for (i = 0; i < N; i++)
    {
        fscanf(fin, "%d", &idx);
        fscanf(fin, "%lf %lf %lf %lf %d",
               &u[idx], &l[idx], &rhs[idx], &d[idx], &p[idx]);
    }
    fclose(fin);

    // Hines algorithm
    HinesAlgo();

    // report execution time
    printf("Execution time of parallel Hines: %lf s\n", t4 - t1);
    printf("setup time: %lf s\n", t2 - t1);
    printf("backward sweep time: %lf s\n", t3 - t2);
    printf("forward sweep time: %lf s\n", t4 - t3);

    // write result to file
    fout = fopen(argv[2], "w");
    for (i = 0; i < N; i++)
        fprintf(fout, "%d %.8e %.8e %.8e %.8e\n",
                i, u[i], l[i], rhs[i], d[i]);
    fclose(fout);

    delete[] u;
    delete[] l;
    delete[] d;
    delete[] rhs;
    delete[] p;
    return 0;
}
