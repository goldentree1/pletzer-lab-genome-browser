import type { ConfigBuilderOpts } from './types';

const myConfig: { [key: string]: ConfigBuilderOpts } = {
  'P.aeruginosa PAO1': {
    ncbiName: 'GCF_000006765.1',
    firstRegion: 'NC_002516.2',
    data: {
      refSeq: 'refseq.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'S.aureus HG001 (NCTC 8325)': {
    ncbiName: 'GCF_000013425.1',
    firstRegion: 'NC_007795.1',
    data: {
      refSeq: 'refseq.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gz',
      coverage: [],
    },
  },
  'S.aureus USA300LAC': {
    ncbiName: 'GCF_000013465.1',
    firstRegion: 'NC_007793.1',
    data: {
      refSeq: 'refseq.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'P.aeruginosa PA14': {
    ncbiName: 'GCF_000014625',
    firstRegion: 'NC_008463.1',
    data: {
      refSeq: 'GCF_000014625.1_ASM1462v1_genomic.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      // coverage: ['PA14_I_1.bw', 'PA14_Un_1.bw'],
      coverage: [
        [
          'PA14_I_Average.bw',
          'PA14_I_1.reheadered.bw',
          'PA14_I_2.reheadered.bw',
          'PA14_I_3.reheadered.bw',
        ],
        [
          'PA14_Un_Average.bw',
          'PA14_Un_1.reheadered.bw',
          'PA14_Un_2.reheadered.bw',
          'PA14_Un_3.reheadered.bw',
        ],
      ],
    },
    extras: [
      // {
      //   type: 'MultiQuantitativeTrack',
      //   trackId: 'TESTmultiwig-coverage',
      //   name: 'TESTCoverage',
      //   assemblyNames: ['asm'],
      //   category: ['Coverage'],
      //   adapter: {
      //     type: 'MultiWiggleAdapter',
      //     bigWigs: [
      //       '/data/GCF_000014625/PA14_I_1.reheadered.bw',
      //       '/data/GCF_000014625/PA14_I_2.reheadered.bw',
      //       '/data/GCF_000014625/PA14_I_3.reheadered.bw',
      //       '/data/GCF_000014625/PA14_Un_1.reheadered.bw',
      //       '/data/GCF_000014625/PA14_Un_2.reheadered.bw',
      //       '/data/GCF_000014625/PA14_Un_3.reheadered.bw',
      //     ],
      //   },
      //   displays: [
      //     {
      //       type: 'MultiLinearWiggleDisplay',
      //       displayId: 'TESTCoverage_multiwiggle-MultiLinearWiggleDisplay',
      //       renderer: {
      //         type: 'XYPlotRenderer',
      //       },
      //       scaleType: 'log',
      //       autoscale: 'global',
      //     },
      //   ],
      // },
    ],
  },
  'K.pneumoniae KPNIH1': {
    ncbiName: 'GCF_000281535.2',
    firstRegion: 'NZ_CP008827.1',
    data: {
      refSeq: 'refseq.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'A.baumannii AB5075-UW': {
    ncbiName: 'GCF_031932345.1',
    firstRegion: 'NZ_JAVSCP010000001.1',
    data: {
      refSeq: 'refseq.fna.gz',
      genomic: 'genomic.gff.sorted.noregion.gff.gz',
      coverage: [],
    },
  },
  'P.aeruginosa LESB58': {
    ncbiName: 'GCF_000026645.1',
    firstRegion: 'NC_011770.1',
    data: {
      refSeq: 'LESB58_ASM2664v1.fasta.gz',
      genomic: 'LESB58_ASM2664v1.gff.sorted.gff.noregion.gz',
      coverage: [['bigwig-bin5.bw'], ['bigwig-bin100.bw']],
    },
    // extras: [
    //   // temporary - extras may be useful
    //   {
    //     type: 'VariantTrack',
    //     trackId: 'VcfVariantTrack1',
    //     name: 'Variants',
    //     assemblyNames: ['asm'],
    //     adapter: {
    //       type: 'VcfTabixAdapter',
    //       vcfGzLocation: {
    //         uri: '/data/GCF_000026645.1/outputCorrected.vcf.gz',
    //         locationType: 'UriLocation',
    //       },
    //       index: {
    //         location: {
    //           uri: '/data/GCF_000026645.1/outputCorrected.vcf.gz.tbi',
    //           locationType: 'UriLocation',
    //         },
    //       },
    //     },
    //   },
    //   {
    //     type: 'QuantitativeTrack',
    //     trackId: 'BigWig-bin5',
    //     name: 'Coverage bin5',
    //     assemblyNames: ['asm'],
    //     adapter: {
    //       type: 'BigWigAdapter',
    //       bigWigLocation: {
    //         uri: '/data/GCF_000026645.1/bigwig-bin5.bw',
    //         locationType: 'UriLocation',
    //       },
    //     },
    //     displays: [
    //       {
    //         type: 'LinearWiggleDisplay',
    //         displayId: 'BigWig-bin5-LinearWiggleDisplay',
    //         height: 150,
    //         renderer: {
    //           type: 'XYPlotRenderer',
    //           color: '#1f77b4',
    //           lineWidth: 1.5,
    //         },
    //         scaleType: 'linear',
    //         autoscale: 'local',
    //       },
    //     ],
    //   },
    //   {
    //     type: 'QuantitativeTrack',
    //     trackId: 'BigWig-bin100',
    //     name: 'Coverage bin100',
    //     assemblyNames: ['asm'],
    //     adapter: {
    //       type: 'BigWigAdapter',
    //       bigWigLocation: {
    //         uri: '/data/GCF_000026645.1/bigwig-bin100.bw',
    //         locationType: 'UriLocation',
    //       },
    //     },
    //     displays: [
    //       {
    //         type: 'LinearWiggleDisplay',
    //         displayId: 'BigWig-bin100-LinearWiggleDisplay',
    //         height: 150,
    //         renderer: {
    //           type: 'XYPlotRenderer',
    //           color: '#d62728',
    //           lineWidth: 2,
    //         },
    //         scaleType: 'linear',
    //         autoscale: 'local',
    //       },
    //     ],
    //   },
    // ],
  },
};

export default myConfig;
