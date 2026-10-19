# 과제 제출 경로와 Git 이력

- 기존 저장소: https://github.com/Gyu-kite/Land_Atmosphere_Group
- 제출 브랜치: `assignment-mmda-20261019`
- 제출 폴더: `MMDA/`
- 폴더 링크: https://github.com/Gyu-kite/Land_Atmosphere_Group/tree/assignment-mmda-20261019/MMDA
- 보고서: [05주차_과제_보고서.md](docs/assignement/05주차_과제_보고서.md)

## 수정 전후

1. `0244d65`: `feat: add audit-da-comparison skill under MMDA`
   - 수정 전 Skill과 research reference 두 파일.
2. 후속 `feat: publish MMDA assignment evidence and documentation` commit:
   - 출력 순서 한 줄 수정, 보고서·A/B/C·재시험 근거·README 사본.

현재 연구 저장소의 최초 commit `c99e648`과 공개 브랜치의 `0244d65`는
조상 및 파일 경로가 다르지만 수정 전 Skill 내용은 동일하다.
관련 없는 연구 commit을 push하지 않기 위해 별도 이력을 만들었다.

```bash
git log --oneline -- MMDA/.agents/skills/audit-da-comparison/SKILL.md
git diff 0244d65 HEAD -- MMDA/.agents/skills/audit-da-comparison/SKILL.md
```

## 범위와 검토

기존 main에서 추가되는 경로는 MMDA/ 아래뿐이다. 원본 연구 파일은 삭제하거나
이동하지 않았다. README는 원본과 동일한 사본이며 내부 상대 경로 일부는
원래 연구 폴더 기준이다. local notebook/cache 없이는 전체 수치 재현은 불가능하다.
과제 원본 PDF, 대용량 모델 output, 위성/ISMN 원자료는 업로드하지 않는다.

최종 파일 검토는 AI가 수행했다. 사용자의 직접 검토·편집 완료로 바꾸어
기록하지 않았으며, 새 대화 독립성 및 한 줄 수정의 일반적 인과 효과도
입증되지 않은 한계로 보고서에 남겼다.

원격 main 자체를 갱신하지 않고 제출 브랜치만 push한다. 따라서 GitHub에서
기본 main 대신 위 제출 브랜치 링크를 사용한다.
