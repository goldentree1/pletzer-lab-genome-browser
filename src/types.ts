import type { createViewState } from '@jbrowse/react-linear-genome-view2';

export type JBrowseConfig = Parameters<typeof createViewState>[0];
export type MyJBrowseConfig = Omit<JBrowseConfig, 'assembly' | 'tracks'>;
export type ViewModel = ReturnType<typeof createViewState>;

// This is for if I make my own Jbrowse conf auto-creator thingy.
export interface JBrowseConfigBuilderOpts {
  assemblyName: string;
  dataDirName: string;

  fastaRefSeqUri: string;
  gff3GenesUri: string;
  vcfVariantsUri?: string;
  bwGeneCountsUri?: string;
  bamUri?: string;
}
