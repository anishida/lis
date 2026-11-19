#!/usr/bin/env python3
"""Property test the actual bounded top-k functions against full qsort.

Extracts the helper functions from the source file under test rather than
copying a second implementation. Runs deterministic tie-rich cases with
both 32-bit and 64-bit LIS_INT models; supports optional ASan/UBSan.
"""
import argparse
from pathlib import Path
import os
import subprocess
import sys
import tempfile

C_TEST = r'''
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef TEST_LIS_INT64
typedef int64_t LIS_INT;
#else
typedef int32_t LIS_INT;
#endif
typedef double LIS_REAL;
typedef double LIS_SCALAR;

/* Exact source selection helpers from src/precon/lis_precon_ilutp.inc: */
@@TOPK@@

static uint32_t state = 0x9e3779b9u;
static uint32_t next_random(void)
{
    uint32_t x = state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    state = x;
    return x;
}
#define LIMIT 400
static LIS_ILUTP_ENTRY input[LIMIT], selected[LIMIT], reference[LIMIT];
static unsigned long long tested = 0;

static int check_case(int m, int cap, int pattern)
{
    int i, kept;
    for(i=0; i<m; i++)
    {
        /* Distinct columns and deliberate ties; preserve unique values
         * even for equal magnitudes, to catch wrong row/column selection. */
        uint32_t r = next_random();
        double v;
        switch(pattern % 6)
        {
        case 0: v = (double)(r % 3); break;
        case 1: v = (double)(r % 17) / 4.; break;
        case 2: v = (double)(i % 5); break;
        case 3: v = (double)(r % 1000) * 1e-12; break;
        case 4: v = ((r & 3u) == 0) ? 0.0 : (double)(r % 101); break;
        default: v = (double)(r % 23) / 13.; break;
        }
        input[i].col = (LIS_INT)i;
        input[i].magnitude = v;
        input[i].value = (r & 0x4000u) ? -v : v;
    }
    for(i=m-1;i>0;i--)
    {
        int j=(int)(next_random()%(uint32_t)(i+1));
        LIS_ILUTP_ENTRY tmp=input[i];
        input[i]=input[j];
        input[j]=tmp;
    }
    memcpy(selected,input,(size_t)m*sizeof(input[0]));
    memcpy(reference,input,(size_t)m*sizeof(input[0]));
    kept=(int)lis_ilutp_keep_largest(selected,(LIS_INT)m,(LIS_INT)cap);
    qsort(reference,(size_t)m,sizeof(reference[0]),lis_ilutp_compare_desc);
    if(kept!=(m<cap?m:cap))
    {
        fprintf(stderr,"FAIL count m=%d cap=%d kept=%d\n",m,cap,kept);
        return 1;
    }
    for(i=0;i<kept;i++)
    {
        if(selected[i].col!=reference[i].col ||
           memcmp(&selected[i].value,&reference[i].value,sizeof(LIS_SCALAR)) ||
           memcmp(&selected[i].magnitude,&reference[i].magnitude,sizeof(LIS_REAL)))
        {
            fprintf(stderr,"FAIL item m=%d cap=%d i=%d (%lld vs %lld)\n",
                    m,cap,i,(long long)selected[i].col,(long long)reference[i].col);
            return 1;
        }
    }
    tested++;
    return 0;
}
int main(void)
{
    int m,cap,k;
    const int sizes[]={0,1,2,3,4,7,8,9,15,16,17,31,32,33,63,64,65,127,128,129,255,256,399,400};
    for(k=0;k<6;k++)
        for(size_t s=0;s<sizeof(sizes)/sizeof(sizes[0]);s++)
        {
            m=sizes[s];
            for(cap=0;cap<=m+2;cap++)
                if(check_case(m,cap,k))return 1;
        }
    for(k=0;k<90000;k++)
    {
        m=(int)(next_random()%(LIMIT+1));
        cap=(int)(next_random()%(LIMIT+1));
        if(check_case(m,cap,k))return 1;
    }
    printf("PASS: %llu deterministic top-k/qsort equivalence cases (%s)\n",
           tested,
#ifdef TEST_LIS_INT64
           "64-bit LIS_INT"
#else
           "32-bit LIS_INT"
#endif
    );
    return 0;
}
'''


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=Path(__file__).resolve().parent.parent / 'src/precon/lis_precon_ilutp.inc')
    p.add_argument('--cc', default=os.environ.get('CC','cc'))
    p.add_argument('--sanitize', action='store_true', help='compile the isolated comparator with ASan and UBSan')
    args = p.parse_args()
    source = args.source.read_text(encoding='utf-8')
    start_marker = 'typedef struct\n{\n    LIS_INT col;'
    end_marker = '/* A min-heap avoids O(lower_fill^2) selection scans.'
    if source.count(start_marker)!=1 or source.count(end_marker)!=1:
        raise SystemExit('ERROR: ILUTP selection source anchors changed; update the property test')
    excerpt = source[source.index(start_marker):source.index(end_marker)]
    if excerpt.count('lis_ilutp_keep_largest(')!=1 or excerpt.count('lis_ilutp_worst_heap_down(')<2:
        raise SystemExit('ERROR: selected helper functions are missing')
    with tempfile.TemporaryDirectory(prefix='ilutp-topk-check-') as folder:
        folder=Path(folder)
        test_src=folder/'selection_test.c'
        test_src.write_text(C_TEST.replace('@@TOPK@@',excerpt),encoding='utf-8')
        for large_int in (False,True):
            exe=folder/('selection_64' if large_int else 'selection_32')
            cmd=[args.cc,'-std=c99','-Wall','-Wextra','-Werror','-O1' if args.sanitize else '-O2','-g',str(test_src),'-o',str(exe)]
            if large_int:cmd+=['-DTEST_LIS_INT64=1']
            if args.sanitize:cmd+=['-fno-omit-frame-pointer','-fsanitize=address,undefined']
            subprocess.run(cmd,check=True)
            env=os.environ.copy()
            if args.sanitize:env['ASAN_OPTIONS']='detect_leaks=1'
            subprocess.run([str(exe)],check=True,env=env)
    return 0

if __name__=='__main__':
    sys.exit(main())
