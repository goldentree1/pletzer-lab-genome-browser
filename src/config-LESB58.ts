import type { JBrowseConfig } from './types';

const config: JBrowseConfig = {
  assembly: {
    name: 'LESB58',
    sequence: {
      type: 'ReferenceSequenceTrack',
      trackId: 'RefSeqTrack',
      adapter: {
        type: 'BgzipFastaAdapter',
        uri: '/data/LESB58/LESB58_ASM2664v1.fasta.gz',
      },
    },
  },
  tracks: [
    {
      type: 'FeatureTrack',
      trackId: 'GFF3GeneTrack',
      name: 'Genes',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'Gff3TabixAdapter',
        gffGzLocation: {
          uri: '/data/LESB58/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/LESB58/LESB58_ASM2664v1.gff.sorted.gff.noregion.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
    },
    {
      type: 'VariantTrack',
      trackId: 'VcfVariantTrack1',
      name: 'Variants',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'VcfTabixAdapter',
        vcfGzLocation: {
          uri: '/data/LESB58/outputCorrected.vcf.gz',
          locationType: 'UriLocation',
        },
        index: {
          location: {
            uri: '/data/LESB58/outputCorrected.vcf.gz.tbi',
            locationType: 'UriLocation',
          },
        },
      },
    },
    {
      type: 'QuantitativeTrack',
      trackId: 'BigWigFromBAMTrack1',
      name: 'Coverage',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'BigWigAdapter',
        bigWigLocation: {
          uri: '/data/LESB58/BigWig_From_BAM-binsize5.bw',
          locationType: 'UriLocation',
        },
      },
      displays: [
        {
          type: 'LinearWiggleDisplay',
          displayId: 'BigWigFromBAMTrack1-LinearWiggleDisplay',
          height: 225,
          // scaleType: 'log', // 'linear' | 'log'
          // autoscale: 'local', // 'local' | 'global'
          // summaryScoreMode: 'mean', // mean | max | min | total
          // logScaleBase: 10, // IS THIS legit!?
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
        tracks: [
          'RefSeqTrack',
          'BigWigFromBAMTrack1',
          'GFF3GeneTrack',
          'VcfVariantTrack1',
        ],
      },
    },
  },

  aggregateTextSearchAdapters: [
    {
      type: 'TrixTextSearchAdapter',
      textSearchAdapterId: 'lesb58-index',
      ixFilePath: {
        uri: '/data/LESB58/trix/assembly1.ix',
        locationType: 'UriLocation',
      },
      ixxFilePath: {
        uri: '/data/LESB58/trix/assembly1.ixx',
        locationType: 'UriLocation',
      },
      metaFilePath: {
        uri: '/data/LESB58/trix/assembly1_meta.json',
        locationType: 'UriLocation',
      },
      assemblyNames: ['LESB58'],
    },
  ],
};

export default config;
