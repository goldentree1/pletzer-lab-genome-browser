import type { createViewState } from '@jbrowse/react-linear-genome-view2';

export type JBrowseConfig = Parameters<typeof createViewState>[0];

export type ViewModel = ReturnType<typeof createViewState>;

/** My personal config - much, much smaller.
 * Used for building full config for JBrowse */
export interface ConfigBuilderOpts {
  ncbiName: string; // e.g., 'GCF_000014625'
  dataDir?: string; // default '/data/$ncbiName/'
  firstRegion: string; // e.g., 'NC_008463.1'
  data: {
    refSeq: string; // default 'REFSEQ.faa.gz'
    genomic: string; // default 'GENOME.gff'
    coverage: string[][];
  };
  extras?: JBrowseConfig['tracks'];
}
