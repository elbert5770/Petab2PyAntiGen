from petab2pyantigen.sbml_to_antimony import clean_antimony

RAW = '''// Created by libAntimony v3.1.3
model *M()

  // Reactions:
  v1_v_0: A => B; k*A;

  // Display Names:
  time_unit is "time";
  v1_v_0 is "v_0";

  model isDescribedBy "http://identifiers.org/doi/10.1/x"
  model created "2019-11-29T09:41:48Z"
  model creator1.givenName "Y"

  // Notes:
  model notes ```
text
more text
```
end

M is "M"
'''


def test_display_names_annotations_and_notes_are_removed():
    text, removed, leftover = clean_antimony(RAW)
    assert " is " not in text and "notes" not in text and "creator1" not in text
    assert "v1_v_0: A => B; k*A;" in text
    assert removed == {"display names": 3, "annotations": 3, "notes": 4}
    assert leftover == []


def test_identifier_is_is_reported_not_rewritten():
    text, _, leftover = clean_antimony("model *T()\n  m1.A is A;\nend\n")
    assert leftover == ["m1.A is A;"]
    assert "m1.A is A" in text
