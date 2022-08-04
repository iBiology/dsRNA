#!/usr/bin/env python

"""
f2b: a pipeline for processing FASTQ files from RNA-Seq to get genomic mapped reads.
"""

import os
import argparse
import json

import cmder
import pandas as pd
from seqflow import Flow, task, logger

from bokeh.models import ColumnDataSource, FactorRange, NumeralTickFormatter, Legend
from bokeh.plotting import figure
from bokeh.models.widgets import Tabs, Panel
from bokeh.embed import components

from jinja2 import Environment, FileSystemLoader, select_autoescape


parser = argparse.ArgumentParser(description=__doc__, prog='dsrna-te')
parser.add_argument('fastq', nargs='+', type=cmder.filename, help='Path to FASTQ file(s), required.')
parser.add_argument('-s', '--sample_sheet', type=cmder.filename, help='Path to sample sheet CSV file, required.')
parser.add_argument('-g', '--genome', type=cmder.dirname,
                    help='Path to STAR genome index directory, optional, default: %(default)s.',
                    default='/project/xchen/genome/mouse/star.index')
parser.add_argument('-t', '--genome_gtf', type=cmder.filename,
                    help='Path to genome GTF file, optional, default: %(default)s.',
                    default='/project/xchen/genome/mouse/GRCm38.p5.gencode.vM16.primary_assembly.annotation.gtf')
parser.add_argument('-T', '--te_gtf', type=cmder.filename,
                    help='Path to TE GTF file, optional, default: %(default)s.',
                    default='/project/xchen/genome/mouse/GRCm38_GENCODE_rmsk_TE.gtf')
parser.add_argument('-o', '--outdir', help='Path to output directory, optional, default: current working directory.')
parser.add_argument('-n', '--cpus', help='Number of CPUs cores, default: %(default)s.', type=int, default=8)
parser.add_argument('-d', '--dryrun', action='store_true',
                    help='Only print out analytical steps with actually processing the data.')
parser.add_argument('-D', '--debug', action='store_true',
                    help='Invoke debug mode that print out more messages and keep intermedia files.')

args = parser.parse_args()
setattr(args, 'outdir', args.outdir or os.getcwd())
setattr(args, 'genome', '/project/xchen/genome/mouse/chr19/star.index' if args.debug else args.genome)
os.makedirs(args.outdir, exist_ok=True)
os.chdir(args.outdir)

if args.sample_sheet:
    SS = pd.read_csv(args.sample_sheet, comment='#')
    ID2NAME = {row.uID: row.name for row in SS.itertuples()}
    FASTQ2UID = {fq: ID2NAME[os.path.basename(os.path.dirname(fq))]
                 for fq in args.fastq if os.path.basename(os.path.dirname(fq)) in ID2NAME}
else:
    ID2NAME = {fq: os.path.basename(fq).split('_')[0] for fq in args.fastq}
    FASTQ2UID = {fq: ID2NAME[fq] for fq in args.fastq}
    
FASTQS, IDS = list(FASTQ2UID.keys()), list(FASTQ2UID.values())

XML = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="mm10" hasGeneTrack="true" hasSequenceTrack="true" version="8">
    <Resources>
        {resources}
    </Resources>
    <Panel height="494" name="DataPanel" width="1263">
        {tracks}
    </Panel>
    <Panel height="607" name="FeaturePanel" width="1263">
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" height="60"
        id="Reference sequence" name="Reference sequence" visible="true"/>
        <Track altColor="0,0,178" attributeKey="Refseq genes" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178"
        colorScale="ContinuousColorScale;0.0;406.0;255,255,255;0,0,178" fontSize="10" height="53" id="mm10_genes"
        name="Refseq genes" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.44765342960288806"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>
"""
RESOURCE = '<Resource path="{bw}"/>'
TRACK = '''
<Track attributeKey="{bw}" autoScale="false" autoscaleGroup="3" fontSize="14" height="60"
    id="{bw}" name="{name}" renderer="BAR_CHART"
    visible="true" windowFunction="mean">
</Track>
'''


@task(inputs=FASTQS, outputs=lambda i: f'{FASTQ2UID[i]}.r1.fastq.gz')
def qc_and_cut_adapter(fastq, out):
    uid = out.replace('.r1.fastq.gz', '')
    cmd = ['fastp', '-i', fastq, '-o', out,
           '-I', fastq.replace('_1.fq.gz', '_2.fq.gz'), '-O', f'{uid}.r2.fastq.gz',
           '--detect_adapter_for_pe', '--dedup', '--thread', args.cpus,
           '-h', f'{uid}.fastp.html', '-j', f'{uid}.fastp.json', '2>', '/dev/null']
    cmder.run(cmd, msg=f'Running QC and adapters removal for {fastq} ...',
              exit_on_error=True, debug=args.debug, pmt=args.debug)
    
    
@task(inputs=qc_and_cut_adapter, outputs=lambda i: i.replace('.r1.fastq.gz', '.unique.r1.fastq.gz'))
def deduplicate_fastq(fastq, out):
    r2, out2 = fastq.replace('.r1.', '.r2.'), out.replace('.r1.', '.r2.')
    cmd = f'fastuniq -i {fastq} {r2} -o {out} -p {out2}'
    cmder.run(cmd, msg=f'Running duplicates from {fastq} ...', exit_on_error=True, debug=args.debug, pmt=args.debug)


@task(inputs=deduplicate_fastq, outputs=lambda i: i.replace('.unique.r1.fastq.gz', '.te.bam'))
def map_reads_to_te(fastq, bam):
    prefix = f'{fastq.replace(".unique.r1.fastq.gz", ".star.te.map")}'
    n = sum([int(os.path.exists(f'{uid}.te.bam')) for uid in IDS])
    genome_load = 'LoadAndRemove' if n == len(args.fastq) - 1 else 'LoadAndKeep'
    try:
        cmder.run(f'[ -d {prefix} ] || mkdir {prefix}')
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', args.cpus,
               '--alignEndsType', 'EndToEnd',
               '--genomeDir', args.genome,
               '--genomeLoad', genome_load,
               '--outBAMcompression', 10,
               '--outFileNamePrefix', f"{prefix}/",
               '--outFilterMultimapNmax', 100,
               '--outFilterMultimapScoreRange', 1,
               '--outSAMmultNmax', 1,
               '--outFilterScoreMin', 10,
               '--outFilterType', 'BySJout',
               '--outReadsUnmapped', 'Fastx',
               '--outSAMattrRGline', 'ID:foo',
               '--outSAMattributes', 'All',
               '--outSAMmode', 'Full',
               '--outSAMtype', 'BAM', 'Unsorted',
               '--outSAMunmapped', 'None',
               '--outStd', 'Log',
               '--readFilesCommand', 'zcat',
               '--readMapNumber', 1_000_000 if args.debug else -1,
               '--readFilesIn', fastq, fastq.replace('.r1.fastq.gz', '.r2.fastq.gz'),
               '>', '/dev/null']
        cmder.run(cmd, msg=f'Mapping reads in {fastq} to TE ...', exit_on_error=True, pmt=args.debug, debug=args.debug)
        cmder.run(f'mv {prefix}/Log.final.out {prefix}.log')
        cmder.run(f'samtools sort -@ {args.cpus} -m 2G -o {bam} {prefix}/Aligned.out.bam', debug=args.debug)
        cmder.run(f'samtools index {bam}')
        cmder.run(f'mv {prefix}/Unmapped.out.mate1 {bam[:-7]}.te.unmap.r1.fastq')
        cmder.run(f'mv {prefix}/Unmapped.out.mate2 {bam[:-7]}.te.unmap.r2.fastq')
    finally:
        if args.debug:
            pass
        else:
            cmder.run(f'rm -r {prefix}')


@task(inputs=[], parent=deduplicate_fastq, outputs=[f'{uid}.bam' for uid in IDS])
def map_reads_to_genome(fastq, bam):
    prefix = f'{bam[:-4]}.genome.map'
    n = sum([int(os.path.exists(f'{name}.bam')) for name in IDS])
    genome_load = 'LoadAndRemove' if n == len(FASTQS) - 1 else 'LoadAndKeep'
    try:
        cmder.run(f'[ -d {prefix} ] || mkdir {prefix}')
        cmd = ['STAR',
               '--runMode', 'alignReads',
               '--runThreadN', args.cpus,
               '--alignEndsType', 'EndToEnd',
               '--genomeDir', args.genome,
               '--genomeLoad', genome_load,
               '--outBAMcompression', 10,
               '--outFileNamePrefix', f"{prefix}/",
               '--outFilterMultimapNmax', 1,
               '--outFilterMultimapScoreRange', 1,
               '--outFilterScoreMin', 10,
               '--outFilterType', 'BySJout',
               '--outReadsUnmapped', 'Fastx',
               '--outSAMattrRGline', 'ID:foo',
               '--outSAMattributes', 'All',
               '--outSAMmode', 'Full',
               '--outSAMtype', 'BAM', 'Unsorted',
               '--outSAMunmapped', 'None',
               '--outStd', 'Log',
               '--readFilesIn', f'{bam[:-4]}.te.unmap.r1.fastq', f'{bam[:-4]}.te.unmap.r2.fastq']
        
        cmder.run(cmd, msg=f'Map repeat elements unmapped reads in {mate} to reference genome ...',
                  pmt=True, debug=args.debug)
        cmder.run(f'mv {prefix}/Log.final.out {bam.replace(".genome.map.bam", ".genome.map.log")}')
        cmder.run(f'samtools sort -@ {args.cpus} -m 2G -o {bam} {prefix}/Aligned.out.bam', debug=args.debug)
        cmder.run(f'samtools index {bam}')
    finally:
        if args.debug or args.keep:
            pass
        else:
            cmder.run(f'rm -r {prefix}')
            
            
@task(inputs=[], outputs=['te.feature.count.csv'], parent=map_reads_to_te)
def te_feature_count(inputs, outputs):
    bam, out = ' '.join([f'{uid}.te.bam' for uid in IDS]), 'featureCounts'
    cmd = f'featureCounts -M -p --countReadPairs -s 2 -T -a {args.te_gtf} {args.cpus} -o {out} {bam}'
    cmder.run(cmd, msg=f'Counting TE features ...', exit_on_error=True, debug=args.debug, pmt=args.debug)
    df = pd.read_csv(out, usecols=lambda c: c == 'Geneid' or c.endswith('bam'), comment='#', dtype=str)
    df = df.rename(columns={c: os.path.basename(c).replace('.bam', '') for c in df.columns})
    df.to_csv(outputs, index=False)
    cmder.run(f'mv {out} te.feature.count.tsv')
    cmder.run(f'mv {out}.summary te.feature.count.summary')


@task(inputs=[], outputs=['genome.feature.count.csv'], parent=map_reads_to_genome)
def genome_feature_count(inputs, outputs):
    bam, out = ' '.join([f'{uid}.bam' for uid in IDS]), 'featureCounts'
    cmd = f'featureCounts -p --countReadPairs --primary -s 2 -T -a {args.genome_gtf} {args.cpus} -o {out} {bam}'
    cmder.run(cmd, msg=f'Counting features ...', exit_on_error=True, debug=args.debug, pmt=args.debug)
    df = pd.read_csv(out, usecols=lambda c: c == 'Geneid' or c.endswith('bam'), comment='#', dtype=str)
    df = df.rename(columns={c: os.path.basename(c).replace('.bam', '') for c in df.columns})
    df.to_csv(outputs, index=False)
    cmder.run(f'mv {out} genome.feature.count.tsv')
    cmder.run(f'mv {out}.summary genome.feature.count.summary')


@task(inputs=map_reads_to_genome, outputs=lambda i: i.replace('.bam', '.bw'))
def bam_to_bw(bam, bw):
    cmd = f'bamCoverage --numberOfProcessors {args.cpus} --effectiveGenomeSize 2652783500 ' \
          f'-b {bam} --normalizeUsing CPM -o {bw}'
    cmder.run(cmd, msg=f'Making {bw} using {bam} ...')


@task(inputs=[], parent=bam_to_bw, outputs=['session.xml'])
def session(inputs, outputs):
    resources, tracks = [], []
    for uid in IDS:
        bam, bw = f'{uid}.bam', f'{uid}.bam'
        resources.append(RESOURCE.format(bw=bw))
        tracks.append(TRACK.format(bw=bw, name=uid.replace('_', ' ')))
    
    with open(outputs, 'w') as o:
        o.write(XML.format(resources='\n'.join(resources), tracks='\n'.join(tracks)))
        
      
# @task(inputs=[], parent=align_reads_to_genome, outputs=['qc.summary.csv'])
def qc_summary(inputs, outputs):
    def fastp_summary(log):
        sample = os.path.basename(log).replace('.fastp.json', '')
        with open(log) as f:
            data = json.load(f)
            result = data['filtering_result']
            reads = sum(result.values())
            result['uID'] = sample
            duplicate = int(reads * data['duplication']['rate'])
            duplication = {'uID': sample, 'unique': reads - duplicate, 'duplication': duplicate}
        return result, duplication
            
    def star_map_count(log):
        counts = [log.replace('.star.genome.map.log', '')]
        with open(log) as f:
            for line in f:
                if 'Number of input reads' in line:
                    reads = int(line.strip().split('\t')[-1])
                elif 'Uniquely mapped reads number' in line:
                    counts.append(int(line.strip().split('\t')[-1]))
                elif 'Number of reads mapped to multiple loci' in line:
                    counts.append(int(line.strip().split('\t')[-1]))
                elif 'Number of reads mapped to too many loci' in line:
                    counts.append(int(line.strip().split('\t')[-1]))
                elif '% of reads unmapped: too many mismatches' in line:
                    counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
                elif '% of reads unmapped: too short' in line:
                    counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
                elif '% of reads unmapped: other' in line:
                    counts.append(int(float(line.strip().split('\t')[-1].replace('%', '')) * reads / 100))
        return counts

    def count_plot(counts):
        def count_bar_plot(data, percent=False):
            samples, categories = data['uID'], [c for c in data.columns if c != 'uID']
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
                      '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
            colors = colors[:len(categories)]
            tooltips = [('uID', '@uID')]
            for category in categories:
                if percent:
                    tooltips.append((f'% {category.lower()}', '@{' + category + '}{0.00 %}'))
                else:
                    tooltips.append((f'# {category.lower()}', '@{' + category + '}{0.00 a}'))
            height = 100 if 15 * len(samples) < 100 else 12 * len(samples)
            print(height, len(samples))
            fig = figure(y_range=FactorRange(factors=samples[::-1]), plot_height=height, sizing_mode='scale_width',
                         tooltips=tooltips)
            fig.add_layout(Legend(), 'right')
            fig.hbar_stack(categories, y='uID', height=0.8, color=colors, legend_label=categories,
                           source=ColumnDataSource(data))
            fig.x_range.start = 0
            fig.x_range.range_padding = 0.1
            formatter = NumeralTickFormatter(format='0%') if percent else NumeralTickFormatter(format='0 a')
            fig.xaxis[0].formatter = formatter
            fig.xaxis.axis_label = 'Percent of reads' if percent else 'Number of reads'
            fig.xaxis.axis_line_color = None
            fig.y_range.range_padding = 0.1
            fig.yaxis.axis_line_color = None
            fig.ygrid.grid_line_color = None
            fig.legend.border_line_color = None
            fig.axis.minor_tick_line_color = None
            fig.axis.major_tick_line_color = None
            fig.outline_line_color = None
            return fig
        
        count_bar = count_bar_plot(counts)
        c = counts.copy().set_index('uID')
        c = c.div(c.sum(axis=1), axis=0).reset_index()
        percent_bar = count_bar_plot(c, percent=True)
        
        count_panel, percent_panel = Panel(child=count_bar, title='Count'), Panel(child=percent_bar, title='Percent')
        return components(Tabs(tabs=[count_panel, percent_panel]))

    reads, duplicates, star = [], [], []
    for uid in IDS:
        read, duplicate = fastp_summary(f'{uid}.fastp.json')
        reads.append(read)
        duplicates.append(duplicate)
        star.append(star_map_count(f'{uid}.star.genome.map.log'))

    reads, scripts = pd.DataFrame(reads).drop_duplicates(), []
    duplicates = pd.DataFrame(duplicates).drop_duplicates()
    columns = ['uID', 'Uniquely mapped reads', 'Reads mapped to multiple loci',
               'Reads mapped to too many loci', 'Reads unmapped: too many mismatches',
               'Reads unmapped: too short', 'Reads unmapped: other']
    star = pd.DataFrame(star, columns=columns).drop_duplicates()
    js, read_div = count_plot(reads)
    scripts.append(js)
    js, dup_div = count_plot(duplicates)
    scripts.append(js)
    js, star_div = count_plot(star)
    scripts.append(js)

    data = {'reads': read_div,
            'duplicate': dup_div,
            'star': star_div,
            'scripts': '\n    '.join(scripts)}

    env = Environment(
        loader=FileSystemLoader('/project/xchen/software/dsrna/data/templates'),
        autoescape=select_autoescape()
    )
    template = env.get_template('qc.summary.html')

    with open(outputs.replace('.csv', '.html'), 'w') as o:
        o.write(template.render(**data))

    with open(outputs, 'w') as o:
        o.write('# Reads\n')
        o.write(f'{reads.to_csv(index=False)}\n')
        o.write('# PCR Duplication\n')
        o.write(f'{duplicates.to_csv(index=False)}\n')
        o.write('# STAR mapping\n')
        o.write(f'{star.to_csv(index=False)}\n')


def cleanup():
    logger.info('Cleaning up ...')
    cmder.run(f"rm {' '.join([f'{uid}.*.fastq*' for uid in IDS])} || true")


@logger.catch()
def main():
    try:
        flow = Flow('dsrna.te', description=__doc__.strip())
        flow.run(dry_run=args.dryrun, cpus=args.cpus)
        if not args.debug and not args.dryrun:
            cleanup()
        if not args.dryrun:
            logger.info('Mission accomplished!')
    finally:
        cmd = """for i in $(ipcs -m | grep "$(whoami)" | awk '{ print $2 }'); do ipcrm -m "$i"; done"""
        cmder.run(cmd, executable='/bin/bash', fmt_cmd=False, log_cmd=False)


if __name__ == '__main__':
    main()
