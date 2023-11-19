from threading import Thread

class MyThread(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs=None, *, daemon=None) -> None:   
        super().__init__(group=group, target=target, 
                        name=name, args=args, 
                        kwargs=kwargs, daemon=daemon)  

    def run(self) -> None:
        try:
            if self._target:
                self.result = self._target(*self._args, **self._kwargs) 
        finally:
            del self._target, self._args, self._kwargs