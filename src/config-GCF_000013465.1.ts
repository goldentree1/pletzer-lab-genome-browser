import { gen } from './config-gen';
// import type { JBrowseConfig } from './types';

// const config: JBrowseConfig = {
//   assembly: {
//     name: 'GCF_000013465.1',
//     sequence: {
//       type: 'ReferenceSequenceTrack',
//       trackId: 'RefSeqTrack',
//       adapter: {
//         type: 'BgzipFastaAdapter',
//         uri: '/data/GCF_000013465.1/refseq.fna.gz',
//       },
//     },
//   },
//   tracks: [
//     {
//       type: 'FeatureTrack',
//       trackId: 'GFF3GeneTrack',
//       name: 'GFF3 Track',
//       assemblyNames: ['GCF_000013465.1'],
//       adapter: {
//         type: 'Gff3TabixAdapter',
//         gffGzLocation: {
//           uri: '/data/GCF_000013465.1/genomic.gff.sorted.noregion.gff.gz',
//           locationType: 'UriLocation',
//         },
//         index: {
//           location: {
//             uri: '/data/GCF_000013465.1/genomic.gff.sorted.noregion.gff.gz.tbi',
//             locationType: 'UriLocation',
//           },
//         },
//       },
//     },
//   ],
//   defaultSession: {
//     name: 'this session',
//     margin: 0,
//     view: {
//       id: 'linearGenomeView',
//       type: 'LinearGenomeView',
//       displayedRegions: [
//         // NC_007793.1	2872769	83	80	81
//         // NC_007790.1	3125	2908860	80	81
//         // NC_007791.1	4439	2912123	80	81
//         // NC_007792.1	37136	2916716	80	81
//         // {
//         //   refName: 'NC_007793.1',
//         //   start: 0,
//         //   end: 6537648,
//         //   assemblyName: 'GCF_000013465.1',
//         // },
//       ],
//       init: {
//         assembly: 'GCF_000013465.1',
//         loc: 'NC_007793.1:1..5,000',
//         tracks: ['RefSeqTrack', 'GFF3GeneTrack'],
//       },
//     },
//   },

//   aggregateTextSearchAdapters: [
//     {
//       type: 'TrixTextSearchAdapter',
//       textSearchAdapterId: 'lesb58-index',
//       ixFilePath: {
//         uri: '/data/GCF_000013465.1/trix/GCF_000013465.1.ix',
//         locationType: 'UriLocation',
//       },
//       ixxFilePath: {
//         uri: '/data/GCF_000013465.1/trix/GCF_000013465.1.ixx',
//         locationType: 'UriLocation',
//       },
//       metaFilePath: {
//         uri: '/data/GCF_000013465.1/trix/GCF_000013465.1_meta.json',
//         locationType: 'UriLocation',
//       },
//       assemblyNames: ['GCF_000013465.1'],
//     },
//   ],
// };

const g = gen({
  ncbiName: 'GCF_000013465.1',
  firstRegion: 'NC_007793.1',
  data: {
    refSeq: 'refseq.fna.gz',
    genomic: 'genomic.gff.sorted.noregion.gff.gz',
    coverage: [],
  },
});

export default g;
