// SPDX-License-Identifier: LGPL-3.0-or-later
//! lsq::pool — the rayon pool one least-squares run executes on.

/// Threads for one optimizer run up to [`LM_BASE_WORK`]; the pool doubles at
/// every further 16× of work. Measured on 32 logical CPUs, ms per LM
/// iteration at 8 / 16 / 32 threads (work = points × layers × columns):
///
/// | case                         | work  | 8     | 16    | 32    |
/// |------------------------------|-------|-------|-------|-------|
/// | 16 films, 61 λ               | 1.1k  | 0.14  | 0.19  | 0.24  |
/// | 80 films, 61 λ               | 5.0k  | 1.36  | 1.46  | 1.66  |
/// | 16 graded films, 61 λ        | 8.4k  | 0.60  | 0.54  | 0.54  |
/// | 80 films, 401 λ              | 33k   | 6.34  | 5.81  | 6.10  |
/// | 32 graded films (1196), 61 λ | 73k   | 5.10  | 3.87  | 3.81  |
/// | 80 graded films (2990), 61 λ | 183k  | 17.3  | 22.4  | 13.8  |
/// | 16 films, 1001 λ × 30 angles | 541k  | 63.0  | 49.1  | 45.5  |
/// | 32 graded films, 201 λ × 30  | 7.3M  | 587   | 381   | 313   |
///
/// A small sweep loses to the cost of waking and coordinating threads that
/// get a few points each; a large one is throughput-bound and takes every
/// thread there is.
pub(crate) const LM_THREADS: usize = 8;

/// Work below which [`LM_THREADS`] is the pool: the 16-film case on 601
/// wavelengths (10.8k) sits at the crossover, a tie between 8 and 16.
pub(crate) const LM_BASE_WORK: usize = 8192;

/// The pool for `work`: [`LM_THREADS`], doubled per 16× of work above
/// [`LM_BASE_WORK`] (8k → 16, 131k → 32, 2.1M → 64), before the cap.
pub(crate) fn lm_threads_for(work: usize) -> usize {
    let mut threads = LM_THREADS;
    let mut bound = LM_BASE_WORK;
    while work >= bound && threads < 1 << 16 {
        threads *= 2;
        bound = bound.saturating_mul(16);
    }
    threads
}

/// Run `f` (one optimizer run) on a rayon worker of a small pool.
///
/// Two costs this avoids, both per parallel call and both paid hundreds of
/// times per run. Called from outside a pool, every sweep is injected into
/// the global queue while the caller sleeps on an OS latch; on a worker it
/// is pushed locally. And a big pool on a small sweep wakes threads that
/// cost more to coordinate than their share saves. The pool has
/// [`lm_threads_for`]`(work)` threads, or `NAVETTE_LM_THREADS` rounded up
/// to a power of two, capped at the global
/// pool's size (so `RAYON_NUM_THREADS` still bounds it); the pools are
/// built once and kept. Results do not depend on the thread count: every
/// sweep collects in index order.
///
/// Already on a worker (a caller that parallelizes over runs), `f` runs
/// where it is, so the caller's own parallelism is left alone.
pub(crate) fn in_lm_pool<R: Send>(work: usize, f: impl FnOnce() -> R + Send) -> R {
    use std::sync::OnceLock;
    const N_POOLS: usize = usize::BITS as usize;
    static POOLS: [OnceLock<Option<rayon::ThreadPool>>; N_POOLS] =
        [const { OnceLock::new() }; N_POOLS];

    if rayon::current_thread_index().is_some() {
        return f();
    }
    let global = rayon::current_num_threads();
    let want = std::env::var("NAVETTE_LM_THREADS")
        .ok()
        .and_then(|s| s.parse::<usize>().ok())
        .unwrap_or_else(|| lm_threads_for(work))
        .max(1)
        .next_power_of_two();
    if want >= global {
        return rayon::scope(|_| f());
    }
    let pool = POOLS[want.trailing_zeros() as usize].get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(want)
            .thread_name(move |i| format!("navette-lm{want}-{i}"))
            .build()
            .ok()
    });
    match pool {
        Some(p) => p.install(f),
        None => rayon::scope(|_| f()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lm_pool_doubles_per_sixteenfold_work() {
        assert_eq!(lm_threads_for(0), 8);
        assert_eq!(lm_threads_for(8191), 8);
        assert_eq!(lm_threads_for(8192), 16);
        assert_eq!(lm_threads_for(131_071), 16);
        assert_eq!(lm_threads_for(131_072), 32);
        assert_eq!(lm_threads_for(2_097_152), 64);
        assert!(lm_threads_for(usize::MAX) <= 1 << 16);
    }
}
