import type { JBrowseCustomConfigHybrid } from './types';
import configJson from '../config.json';

const myConfig: { [key: string]: JBrowseCustomConfigHybrid } = configJson;
export default myConfig;
