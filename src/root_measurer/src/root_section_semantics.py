"""Whole connected-area equivalent radius, distinct from a single-branch radius."""
import math


def connected_section(row):
    area=row.get('area');status=row.get('status','NOT_MEASURED')
    qualified=(status.startswith('CLOSED_DISK_') or status in
               ['MISSING_MIXED_BRANCH_SECTION','NONDISK_SECTION_HOLES_NOT_FILLED'])
    available=qualified and area is not None and math.isfinite(area) and area>0
    label=('mixed_branch' if status=='MISSING_MIXED_BRANCH_SECTION' else
           'section_with_holes' if status=='NONDISK_SECTION_HOLES_NOT_FILLED' else
           'ordinary_section' if available else 'unqualified_section')
    return dict(radius=math.sqrt(area/math.pi) if available else None,
                area=area if available else None,available=available,classification=label,
                holes_filled=False,definition='sqrt(area of entire original connected material section/pi); holes deducted')
