import { useEffect } from 'react';

export function useShiftScrollPan(containerId?: string) {
  useEffect(() => {
    const targetElement: HTMLElement | null = containerId
      ? document.getElementById(containerId)
      : (document.scrollingElement as HTMLElement | null);

    if (!targetElement) return;

    const handler = (e: WheelEvent): void => {
      if (!e.shiftKey) return;

      e.preventDefault();

      // allow both axes while holding shift
      targetElement.scrollLeft += e.deltaY;
      targetElement.scrollTop += e.deltaX;
    };

    targetElement.addEventListener('wheel', handler, { passive: false });

    return () => {
      targetElement.removeEventListener('wheel', handler);
    };
  }, [containerId]);
}
