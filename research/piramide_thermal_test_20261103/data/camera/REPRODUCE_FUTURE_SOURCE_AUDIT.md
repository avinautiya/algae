# Released future-batch source audit

Authorized release follows final prediction publication `ebd47166c8bd1715eff3668337d5f18d4ead4e52`. The receipt is `model/frozen_piramide_future_batch_final/EVALUATION_RELEASE_RECEIPT.json`. The pre-existing role-gated audit script checks the receipt and source hashes before dereferencing authorized endpoints.

From the round2 root, using the isolated parser environment:

```
data/venv312/bin/python data/camera/postscore_source_audit.py data/camera model/frozen_piramide_future_batch_final/EVALUATION_RELEASE_RECEIPT.json --output data/camera/FUTURE_BATCH_ENDPOINT_SOURCE_AUDIT.json
data/venv312/bin/python data/camera/verify_future_unit_audit.py
```

The three source ice-ablation amounts are59.046116588954,61.451154426008 and89.054068967277mm water equivalent. Each independently recomputed raw distance change×917 agrees with supplied ablation change×1000 to floating-point precision. Thus supplied ablation is consistent with metres water equivalent at these released endpoints, despite the README's millimetres label. The frozen raw-distance operator remains unchanged.

Only Dec30 at00:00, Jan13 at03:00, Jan27 at03:00 and Feb10 at03:00 were released. The excluded Feb24 endpoint remains unopened. All three endpoint differences are positive; this does not rule out interior reference resets or quantify manual correction error. No camera-precision bound or confidence interval is established. The source check neither changes the comparison models nor reverses their failed overall-superiority result.
