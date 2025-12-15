import type { createViewState } from '@jbrowse/react-linear-genome-view2';
type JBrowseConfig = Parameters<typeof createViewState>[0];

export const config: JBrowseConfig = {
  assembly: {
    name: 'LESB58',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        // uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz',
        fastaLocation: {
          uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz',
          locationType: 'UriLocation',
        },
        faiLocation: {
          uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz.fai',
          locationType: 'UriLocation',
        },
        gziLocation: {
          uri: '/data/NC011770/LESB58_ASM2664v1.fasta.gz.gzi',
          locationType: 'UriLocation',
        },
      },
    },

    // refNameAliases: {
    //   adapter: {
    //     type: 'RefNameAliasAdapter',
    //     uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt',
    //   },
    // },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'GFF3 Track',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/NC011770/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/NC011770/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
      displays: [
        {
          type: 'LinearBasicDisplay',
          displayId: 'GFF3GeneTrack-LinearBasicDisplay',
          renderer: {
            type: 'SvgFeatureRenderer',
            color1: "jexl:get(feature,'strand')>0?'#00FF00':'red'",
          },
        },
      ],
    },
  ],
  defaultSession: {
    name: 'this session',
    margin: 0,
    view: {
      id: 'linearGenomeView',
      type: 'LinearGenomeView',
      displayedRegions: [
        {
          refName: 'NC_011770.1',
          start: 0,
          end: 6601757,
          assemblyName: 'LESB58',
        },
      ],
      init: {
        assembly: 'LESB58',
        loc: 'NC_011770.1:1..5,000',
        tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
      },
    },
  },

  configuration: {
    rpc: {
      defaultDriver: 'WebWorkerRpcDriver',
    },
  },

  makeWorkerInstance: () => {
    return new Worker(new URL('./rpcWorker', import.meta.url), {
      type: 'module',
    });
  },
};

// export const config: JBrowseConfig = {
//   assembly: {
//     name: 'hg38',
//     sequence: {
//       type: 'ReferenceSequenceTrack',
//       trackId: 'GRCh38-ReferenceSequenceTrack',
//       adapter: {
//         type: 'BgzipFastaAdapter',
//         uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz',
//       },
//     },
//     refNameAliases: {
//       adapter: {
//         type: 'RefNameAliasAdapter',
//         uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/hg38_aliases.txt',
//       },
//     },
//   },
//   tracks: [
//     {
//       type: 'FeatureTrack',
//       trackId: 'genes',
//       name: 'NCBI RefSeq Genes',
//       assemblyNames: ['hg38'],
//       category: ['Genes'],
//       adapter: {
//         type: 'Gff3TabixAdapter',
//         uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/ncbi_refseq/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff.gz',
//       },
//       textSearching: {
//         textSearchAdapter: {
//           type: 'TrixTextSearchAdapter',
//           textSearchAdapterId: 'gff3tabix_genes-index',
//           uri: 'https://jbrowse.org/genomes/GRCh38/ncbi_refseq/trix/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.sorted.gff.gz.ix',
//           assemblyNames: ['GRCh38'],
//         },
//       },
//     },
//     {
//       type: 'FeatureTrack',
//       trackId: 'repeats_hg38',
//       name: 'Repeats',
//       assemblyNames: ['hg38'],
//       category: ['Annotation'],
//       adapter: {
//         type: 'BigBedAdapter',
//         uri: 'https://jbrowse.org/genomes/GRCh38/repeats.bb',
//       },
//     },
//     {
//       type: 'AlignmentsTrack',
//       trackId: 'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome',
//       name: 'NA12878 Exome',
//       assemblyNames: ['hg38'],
//       category: ['1000 Genomes', 'Alignments'],
//       adapter: {
//         type: 'CramAdapter',
//         uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/alignments/NA12878/NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram',

//         sequenceAdapter: {
//           type: 'BgzipFastaAdapter',
//           uri: 'https://jbrowse.org/genomes/GRCh38/fasta/hg38.prefix.fa.gz',
//         },
//       },
//     },
//     {
//       type: 'VariantTrack',
//       trackId:
//         'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf',
//       name: '1000 Genomes Variant Calls',
//       assemblyNames: ['hg38'],
//       category: ['1000 Genomes', 'Variants'],
//       adapter: {
//         type: 'VcfTabixAdapter',
//         uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/variants/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz',
//       },
//     },
//     {
//       type: 'QuantitativeTrack',
//       trackId: 'hg38.100way.phyloP100way',
//       name: 'hg38.100way.phyloP100way',
//       category: ['Conservation'],
//       assemblyNames: ['hg38'],
//       adapter: {
//         type: 'BigWigAdapter',
//         uri: 'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw',
//       },
//     },
//     {
//       type: 'AlignmentsTrack',
//       trackId: 'skbr3_pacbio',
//       name: 'SKBR3 pacbio',
//       assemblyNames: ['hg38'],
//       adapter: {
//         type: 'BamAdapter',
//         uri: 'https://s3.amazonaws.com/jbrowse.org/genomes/GRCh38/skbr3/SKBR3_Feb17_GRCh38.sorted.bam',
//       },
//     },
//   ],
//   defaultSession: {
//     name: 'this session',
//     margin: 0,
//     view: {
//       id: 'linearGenomeView',
//       type: 'LinearGenomeView',
//       init: {
//         assembly: 'hg38',
//         loc: '10:29,838,565..29,838,850',
//         tracks: [
//           'GRCh38-ReferenceSequenceTrack',
//           'NA12878.alt_bwamem_GRCh38DH.20150826.CEU.exome',
//           'hg38.100way.phyloP100way',
//           'ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf',
//         ],
//       },
//     },
//   },
// };
