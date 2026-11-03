# Endpoint repair: taxonomy, not a retroactive scientific win

Pinned actual branch `claude/recovery-m5` commit455bf85856bcf11771f0b37bd63cc58e89f05f32. Four downloaded source files verify against their Git blob hashes; seeSOURCE_MANIFEST.json. Historical scorer bytes are preserved unmodified inhistorical_score_trace.py from the previously supplied packet.

## Diagnosis

The old scorer is internally consistent: its `reconnection` means a rising threshold crossing AFTER the first loss. An initially dry route that wets during a pulse and subsequently dries has its first loss at pulse termination. It therefore receives `lost_without_reconnection`, although it did temporarily wet. This is an insufficiently expressive taxonomy for that scientific question, not evidence that sustained recovery occurred.

The actual frozen protocol predicted NOTsustained for all36evolving cases and finer classes for pulse cases. Primary binary predictions match36/36. Finer label predictions mismatch27/36. Report BOTH. Do not describe the entire frozen primary criterion as failed solely because finer labels mismatch, and do not convert those27 mismatches into classification successes by applyingv2 retrospectively.

An additional actual runner defect is the history identity export: it asks for `preloss_prefix_sha256` or `fingerprint`, neither of which the scorer returns. All72stored IDs are consequentlynull. The actual field is `sampled_preloss_history_sha256`. The new helper fails if it is missing. Equal sampled prefixes only demonstrate duplicate sampled histories; they never prove statistical independence.

## Prospective v2

`score_trace_v2.py` explicitly distinguishes initially disconnected/never wet, initially disconnected/temporarily wet/then dry, initially disconnected/sustained post-pulse access, and initially disconnected/partial post-pulse access. Originally wet histories keep maintained/lost/reconnected/relost classes. It separately records first threshold touch and positive-duration contact. It preserves the old score in every result and never changes a historical prediction rule.

Fourteen regression tests pass, including the reported mismatch fixture, sustained/partial/no-access cases, isolated threshold touch, signed heat, missing history IDs, sampling gaps and nonfinite values. The samepiecewise-linear observation operator is declared explicitly.

**Timing limit:** the ODE runner uses discontinuous pulse forcing but stores one ordinary sample at each switching timestamp. Linear interpolation of the stored response smears a jump over a sample interval. Neither v1 nor thisv2 recovers exact discontinuity timing. For future jump-accurate endpoints, freeze a separate one-sided/event-aware observation contract and retain exact event timestamps or left/right limits. Do not present interpolation-defined crossing times as exact physical pulse-response times. The no-sustained conclusion remains unaffected in these records.

Run `python3 -m unittest -v test_score_trace_v2` and `python3 audit_frozen.py` here. No Claude file, frozen score, prediction or historical result was edited. This is not field validation, novelty certification or a guarantee of positive outcomes.
