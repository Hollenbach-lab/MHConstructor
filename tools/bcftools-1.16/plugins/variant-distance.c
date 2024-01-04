/*
    Copyright (C) 2022 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/
/*
    Add custom INFO field with the distance to the nearest variant

*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <assert.h>
#include "bcftools.h"
#include "filter.h"
#include "vcfbuf.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    char *tag;
    int argc, region_is_file, target_is_file, output_type, clevel;
    int regions_overlap, targets_overlap;
    int32_t prev_dist;
    char **argv, *region, *target, *fname, *output_fname;
    char *filter_str;
    filter_t *filter;
    int filter_logic;    // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    htsFile *out_fh;
    bcf_hdr_t *hdr;
    bcf_srs_t *sr;
    vcfbuf_t *buf;
}
args_t;

const char *about(void)
{
    return "Annotate sites with distance to the nearest variant\n";
}

static const char *usage_text(void)
{
    return
        "\n"
        "About: Annotate sites with distance to the nearest variant\n"
        "Usage: bcftools +variant-distance [Options]\n"
        "Options:\n"
        "   -n, --tag-name STR              The tag name to be added [DIST]\n"
        "Common options:\n"
        "   -e, --exclude EXPR              Exclude sites for which the expression is true\n"
        "   -i, --include EXPR              Include only sites for which the expression is true\n"
        "   -o, --output FILE               Write output to the FILE [standard output]\n"
        "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n"
        "   -r, --regions REGION            Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -t, --targets REGION            Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "Examples:\n"
        "   bcftools +variant-distance input.bcf -Ob -o output.bcf\n"
        "\n";
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->region )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->region, args->region_is_file)<0 ) error("Failed to read the regions: %s\n",args->region);
    }
    if ( args->target && bcf_sr_set_targets(args->sr, args->target, args->target_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->target);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    if ( !args->tag ) args->tag = strdup("DIST");
    bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Integer,Description=\"Distance to the  nearest variant\">",args->tag);

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));

    if ( bcf_hdr_write(args->out_fh, args->hdr)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);

    args->buf = vcfbuf_init(args->hdr, 0);
    vcfbuf_set_opt(args->buf,int,VCFBUF_DUMMY,1)
}
static void destroy_data(args_t *args)
{
    if ( args->filter )
        filter_destroy(args->filter);
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    bcf_sr_destroy(args->sr);
    vcfbuf_destroy(args->buf);
    free(args->tag);
    free(args);
}
static void flush(args_t *args, int n)
{
    int i;
    for (i=0; i<n; i++)
    {
        bcf1_t *rec = vcfbuf_flush(args->buf, 1);
        if ( args->prev_dist )
        {
            if ( bcf_update_info_int32(args->hdr,rec,args->tag,&args->prev_dist,1)!=0 )
                error("[%s] Error: failed to update INFO/%s=%d\n",__func__,args->tag,args->prev_dist);
        }
        if ( bcf_write(args->out_fh,args->hdr,rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    }
}
static void process(args_t *args, int done)
{
    if ( done )
    {
        // If the buffer is not empty, all records must be duplicate positions and
        // the previous site may or may not be set
        flush(args, vcfbuf_nsites(args->buf));
        return;
    }

    // Do filtering when requested
    bcf1_t *rec  = bcf_sr_get_line(args->sr,0);
    if ( args->filter )
    {
        int pass = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_EXCLUDE ) pass = pass ? 0 : 1;
        if ( !pass ) return;
    }
    bcf_sr_t *sr = bcf_sr_get_reader(args->sr,0);
    sr->buffer[0] = vcfbuf_push(args->buf, rec);

    // Find the block of duplicate records (the same chr+position)
    int i, n = vcfbuf_nsites(args->buf);
    bcf1_t *rec0 = vcfbuf_peek(args->buf, 0);
    for (i=1; i<n; i++)
    {
        rec = vcfbuf_peek(args->buf, i);
        if ( rec->rid != rec0->rid ) break;
        if ( rec->pos != rec0->pos ) break;
    }
    if ( i==n ) return;             // cannot flush yet, either a single record or all duplicates
    if ( rec->rid != rec0->rid )    // the new record is a different chr
    {
        flush(args, i);
        args->prev_dist = 0;
        return;
    }

    // the new record is a different pos
    if ( !args->prev_dist || args->prev_dist > rec->pos - rec0->pos ) args->prev_dist = rec->pos - rec0->pos;
    flush(args, i);
    args->prev_dist = rec->pos - rec0->pos;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->output_fname = "-";
    args->clevel = -1;
    static struct option loptions[] =
    {
        {"tag-name",required_argument,NULL,'n'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"targets-overlap",required_argument,NULL,2},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "r:R:t:T:o:O:n:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'n': args->tag = strdup(optarg); break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'T': args->target_is_file = 1; // fall-through
            case 't': args->target = optarg; break;
            case 'R': args->region_is_file = 1; // fall-through
            case 'r': args->region = optarg; break;
            case  2 :
                      args->targets_overlap = parse_overlap_option(optarg);
                      if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                      break;
            case  3 :
                      args->regions_overlap = parse_overlap_option(optarg);
                      if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                      break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      }
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");

    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s",usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s",usage_text());
    else args->fname = argv[optind];

    init_data(args);

    while ( bcf_sr_next_line(args->sr) ) process(args,0);
    process(args,1);

    destroy_data(args);
    return 0;
}


