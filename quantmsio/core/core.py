"""
Core classes for quantmsio provide utilitiy classes for the quantms.io format.
"""

import shutil
import uuid

import diskcache as diskcache


class DiskCache:
    """
    Disk cache class for quantms.io. This class is used to store dictionaries in disks when the size of the dictionary
    is too large to be stored in memory.
    """

    def __init__(self, name_prefix: str):
        # Create a cache name using a hash and uuid
        if name_prefix is None:
            name_prefix = "generic"
        self._cache_name = str("_cache_name_{}_{}".format(name_prefix, uuid.uuid4().hex))
        self.cache = diskcache.Cache(self._cache_name, statistics=True)
        self.cache.create_tag_index()

    def get_item(self, key):
        """
        Get an item from the cache
        """
        return self.cache[key]

    def get_first_subkey(self, subkey):
        """
        Get an item from the cache
        """
        for key in self.cache.iterkeys():
            if subkey in key:
                return self.cache[key]
        return None

    def add_item(self, key, value):
        self.cache[key] = value

    def delete_itm(self, key):
        del self.cache[key]

    def length(self):
        return len(self.cache)

    def get_all_keys(self):
        return self.cache.iterkeys()

    def contains(self, key):
        return self.cache.__contains__(key)

    def close(self):
        """
        Close the cache and delete the cache files. This method should be called when the cache is not needed anymore.
        """
        self.cache.clear()
        self.cache.close()
        shutil.rmtree(self._cache_name)  # remove the disk cache files

    def get_stats(self):
        hits, misses = self.cache.stats()
        count = self.cache._sql("select count(*) from Cache").fetchone()
        return {
            "hits": hits,
            "misses": misses,
            "count": count[0],
            "size": self.cache.volume(),
        }
