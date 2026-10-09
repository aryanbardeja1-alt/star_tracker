// Labels for identified stars.
//
// The catalogue carries an IAU proper name for 335 of its 8763 entries, and
// `MatchDto` now brings it along with the identification, so a named star is
// labelled by name and every other star by its HIP number. The name has to
// travel with the match rather than be looked up here: the tracker lives on a
// worker, so asking per star would mean a round-trip per label.

import type { MatchDto } from "../wasm.js";

/** A short label for an identified star: its proper name, or its HIP number. */
export function starLabel(matched: MatchDto): string {
  return matched.name ?? `HIP ${matched.id}`;
}
