import { useEffect } from "react";

/**
 * Scrolls to hash in url when component mounts and removes hash from url when component unmounts
 * @param valid_hashes - what hashes in url should trigger scroll
 */
const useScrollToHash = (
  valid_hashes: string[],
  removeHash: boolean = true,
) => {
  useEffect(() => {
    const hash = window?.location?.hash.split("#")?.[1];
    if (hash && valid_hashes.includes(hash)) {
      const hashElement = document.getElementById(hash);
      if (hashElement) {
        setTimeout(
          () =>
            hashElement.scrollIntoView({ block: "center", behavior: "smooth" }),
          500,
        );
      }
    }
    // Remove hash when component unmounts
    return () => (removeHash ? history.replaceState(null, "", " ") : undefined);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);
};

export default useScrollToHash;
