import type { createViewState } from '@jbrowse/react-linear-genome-view2';

export type JBrowseConfig = Parameters<typeof createViewState>[0];

export type ViewModel = ReturnType<typeof createViewState>;

export interface JBrowseConfigBuilderOpts {
  assemblyName: string;
  dataDirName: string;

  fastaRefSeqUri: string;
  gff3GenesUri: string;
  vcfVariantsUri?: string;
  bwGeneCountsUri?: string;
  bamUri?: string;
}
