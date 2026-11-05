# 11. The pancake_2layers_180fps.mp4 video

**You still have it. It was renamed, not lost.** The same video is on disk at:

```
archive/code-original/plots/pancake_2layers_n1000_180fps.mp4
```

I compared the git hashes: this file is byte-identical to the original `pancake_2layers_180fps.mp4`
(blob `078312c`, 7.5 MB). I did not play the video. I also did not change, move or restore anything.

## What happened to the file

| Date | Commit | Change |
|---|---|---|
| 2026-06-15 | `8eaea7c` "habemus 2layers" | added as `code/pancake_2layers_180fps.mp4` |
| 2026-06-16 | `18dd7bb` "we are onto something" | renamed to `code/pancake_2layers_n1000_180fps.mp4` |
| 2026-07-24 | `c2ccadc` "move plots" | moved to `code/plots/` |
| 2026-10-06 | `5316b35` "Reorganise repository" | moved to `code-original/plots/` |
| 2026-11-03 | `b21fbb4` "november report" | removed from git tracking |

## Backups

- **The copy in `archive/` is not tracked**, because `.gitignore` excludes both `archive/` and `*.mp4`.
  Deleting that folder would remove this copy.
- **Git history still has the video**, on `dev` and `origin/dev`. To get it back under its original name, run this
  from the repo root:

  ```
  git show 8eaea7c:code/pancake_2layers_180fps.mp4 > pancake_2layers_180fps.mp4
  ```

  Use Git Bash for this. In Windows PowerShell 5.1, `>` can corrupt binary files.

## If this isn't the video you remember

The other videos from that time are in the same folder: `pancake_2layers_n600_180fps.mp4`,
`pancake_3layers_n2000_180fps.mp4` and `pancake_video.mp4`.
