"""Background jobs for the web app: one at a time, since the MLIP holds the CPU."""
import threading
import traceback

from mxbuilder import workflow


class Busy(Exception):
    pass


class Job:
    def __init__(self, job_id, kind):
        self.id, self.kind = job_id, kind
        self.state = "running"          # running | done | stopped | error
        self.events, self.result, self.error = [], None, None
        self.stop = threading.Event()


class Jobs:
    def __init__(self, runs_dir):
        self.runs_dir = runs_dir
        self.jobs = {}
        self.current = None
        self.lock = threading.Lock()

    def start(self, kind, inputs, **options):
        """kind is "run" (workflow.run) or "scan" (workflow.scan, options grid and relax_z)."""
        job_id = inputs.run_id() if kind == "run" else workflow.scan_id(inputs, options["grid"], options["relax_z"])
        with self.lock:
            if self.current and self.current.state == "running":
                raise Busy("Another job is running. Stop it first.")
            job = self.current = self.jobs[job_id] = Job(job_id, kind)
        threading.Thread(target=self._work, args=(job, inputs, options), daemon=True).start()
        return job_id

    def _work(self, job, inputs, options):
        task = workflow.run if job.kind == "run" else workflow.scan
        try:
            job.result = task(inputs, self.runs_dir, progress=job.events.append, should_stop=job.stop.is_set, **options)
            job.state = "stopped" if job.result.get("stopped") else "done"
        except Exception as exc:
            traceback.print_exc()
            job.error, job.state = f"{type(exc).__name__}: {exc}", "error"

    def get(self, job_id):
        return self.jobs.get(job_id)
