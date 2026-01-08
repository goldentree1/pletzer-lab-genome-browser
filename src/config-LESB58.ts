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
      trackId: 'BigWig-bin5',
      name: 'Coverage bin5',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'BigWigAdapter',
        bigWigLocation: {
          uri: '/data/LESB58/bigwig-bin5.bw',
          locationType: 'UriLocation',
        },
      },
      displays: [
        {
          type: 'LinearWiggleDisplay',
          displayId: 'BigWig-bin5-LinearWiggleDisplay',
          height: 150,
          renderer: {
            type: 'XYPlotRenderer',
            color: '#1f77b4',
            lineWidth: 1.5,
          },
          scaleType: 'linear',
          autoscale: 'local',
        },
      ],
    },
    {
      type: 'QuantitativeTrack',
      trackId: 'BigWig-bin100',
      name: 'Coverage bin100',
      assemblyNames: ['LESB58'],
      adapter: {
        type: 'BigWigAdapter',
        bigWigLocation: {
          uri: '/data/LESB58/bigwig-bin100.bw',
          locationType: 'UriLocation',
        },
      },
      displays: [
        {
          type: 'LinearWiggleDisplay',
          displayId: 'BigWig-bin100-LinearWiggleDisplay',
          height: 150,
          renderer: {
            type: 'XYPlotRenderer',
            color: '#d62728',
            lineWidth: 2,
          },
          scaleType: 'linear',
          autoscale: 'local',
        },
      ],
    },
    {
      type: 'MultiQuantitativeTrack',
      trackId: 'Coverage_multiwiggle',
      name: 'Coverage binsize5 vs binsize100',
      assemblyNames: ['LESB58'],
      category: ['Coverage'],
      adapter: {
        type: 'MultiWiggleAdapter',
        bigWigs: [
          '/data/LESB58/bigwig-bin100.bw',
          '/data/LESB58/bigwig-bin5.bw',
        ],
      },
      displays: [
        {
          type: 'MultiLinearWiggleDisplay',
          displayId: 'Coverage_multiwiggle-MultiLinearWiggleDisplay',
          renderer: {
            type: 'XYPlotRenderer',
          },
          scaleType: 'log',
          autoscale: 'global',
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
          'BigWig-bin5',
          'BigWig-bin100',
          'Coverage_multiwiggle',
          // 'BigWig-overlay',
          // 'BigWigFromBAMTrack1',
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
