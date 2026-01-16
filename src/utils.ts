import { useEffect, useState } from 'react';

/** useState<boolean>(), but using localStorage for persistence */
export function useStoredStateBoolean(key: string, defaultValue: boolean) {
  const [value, setValue] = useState<boolean>(() => {
    const stored = localStorage.getItem(key);
    return stored === null ? defaultValue : stored === 'true';
  });

  useEffect(() => {
    localStorage.setItem(key, String(value));
  }, [key, value]);

  return [value, setValue] as const;
}

/** useState<string>(), but using localStorage for persistence */
export function useStoredStateString(key: string, defaultValue: string) {
  const [value, setValue] = useState<string>(() => {
    const stored = localStorage.getItem(key);
    return stored === null ? defaultValue : stored;
  });

  useEffect(() => {
    localStorage.setItem(key, value);
  }, [key, value]);

  return [value, setValue] as const;
}
