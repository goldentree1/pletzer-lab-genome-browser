import type { JBrowseCustomConfig } from './types';
import configJson from './config.json';

const myConfig: { [key: string]: JBrowseCustomConfig } = configJson;
export default myConfig;
