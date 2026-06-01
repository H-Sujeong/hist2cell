#!/usr/bin/env node
/**
 * upload-to-notion.js
 *
 * 로컬 Markdown 파일을 Notion 페이지로 업로드한다.
 *   - md 안의 로컬 이미지(./img/*.png 등)를 Notion File Upload API로 올리고
 *   - @tryfabric/martian 으로 md → Notion 블록 변환
 *   - 이미지 블록을 file_upload 참조로 교체
 *   - 2000자 초과 rich_text 분할
 *   - 부모 페이지 밑에 새 페이지 생성 + children 삽입
 *
 * 사용법:
 *   NOTION_TOKEN=secret_xxx node upload-to-notion.js <md경로> <부모페이지ID>
 *
 * 요구사항: Node 18+ (전역 fetch / FormData / Blob), 그리고
 *   npm i @tryfabric/martian
 */

'use strict';

const fs = require('fs');
const path = require('path');
const { markdownToBlocks } = require('@tryfabric/martian');

// ────────────────────────────────────────────────────────────────────────────
// 설정
// ────────────────────────────────────────────────────────────────────────────
const NOTION_TOKEN = process.env.NOTION_TOKEN;
const NOTION_API = 'https://api.notion.com/v1';

// 일반 블록/페이지 API 버전.
const NOTION_VERSION = '2022-06-28';
// File Upload API 는 비교적 최신에 GA 된 엔드포인트라 별도 상수로 분리해 둔다.
// 현재는 동일 버전("2022-06-28")으로 동작하지만, 추후 Notion 이 버전을 요구하면
// 여기만 바꾸면 된다. (env NOTION_FILE_UPLOAD_VERSION 으로 override 가능)
const FILE_UPLOAD_VERSION =
  process.env.NOTION_FILE_UPLOAD_VERSION || NOTION_VERSION;

const RICH_TEXT_LIMIT = 2000; // Notion rich_text.text.content 최대 길이
const CHILDREN_PER_REQUEST = 100; // children 배열 1회 요청 최대 길이

const IMAGE_EXT_MIME = {
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.jpeg': 'image/jpeg',
  '.gif': 'image/gif',
  '.webp': 'image/webp',
  '.svg': 'image/svg+xml',
  '.bmp': 'image/bmp',
  '.tif': 'image/tiff',
  '.tiff': 'image/tiff',
};

// ────────────────────────────────────────────────────────────────────────────
// 공통 fetch 래퍼
// ────────────────────────────────────────────────────────────────────────────
async function notionFetch(url, { method = 'GET', headers = {}, body, version } = {}) {
  const res = await fetch(url, {
    method,
    headers: {
      Authorization: `Bearer ${NOTION_TOKEN}`,
      'Notion-Version': version || NOTION_VERSION,
      ...headers,
    },
    body,
  });

  const text = await res.text();
  let json;
  try {
    json = text ? JSON.parse(text) : {};
  } catch {
    json = { raw: text };
  }

  if (!res.ok) {
    const msg = json && json.message ? json.message : text;
    throw new Error(`Notion API ${method} ${url} → ${res.status}: ${msg}`);
  }
  return json;
}

// ────────────────────────────────────────────────────────────────────────────
// 1) md 에서 로컬 이미지 경로 추출
//    - ![alt](path)  마크다운 이미지
//    - <img src="path">  HTML 이미지
//    http(s):// 로 시작하면 원격이므로 건너뛴다.
// ────────────────────────────────────────────────────────────────────────────
function extractLocalImagePaths(md) {
  const urls = new Set();

  const mdImg = /!\[[^\]]*\]\(\s*(<[^>]+>|[^)\s]+)(?:\s+["'][^"']*["'])?\s*\)/g;
  const htmlImg = /<img\b[^>]*\bsrc\s*=\s*["']([^"']+)["']/gi;

  let m;
  while ((m = mdImg.exec(md)) !== null) {
    let u = m[1].trim();
    if (u.startsWith('<') && u.endsWith('>')) u = u.slice(1, -1).trim();
    urls.add(u);
  }
  while ((m = htmlImg.exec(md)) !== null) {
    urls.add(m[1].trim());
  }

  return [...urls].filter((u) => u && !/^https?:\/\//i.test(u) && !u.startsWith('data:'));
}

// ────────────────────────────────────────────────────────────────────────────
// 2) 이미지 1개를 Notion File Upload API 로 업로드
//    (a) POST /v1/file_uploads        → { id, upload_url }
//    (b) POST {upload_url} (= .../send) 로 바이너리 multipart 전송
//    반환: file_upload id
// ────────────────────────────────────────────────────────────────────────────
async function uploadImage(absPath) {
  const filename = path.basename(absPath);
  const ext = path.extname(absPath).toLowerCase();
  const mime = IMAGE_EXT_MIME[ext] || 'application/octet-stream';
  const buffer = fs.readFileSync(absPath);

  // (a) 업로드 객체 생성 (단일 파트, <20MB 기준)
  const created = await notionFetch(`${NOTION_API}/file_uploads`, {
    method: 'POST',
    version: FILE_UPLOAD_VERSION,
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      mode: 'single_part',
      filename,
      content_type: mime,
    }),
  });

  const uploadUrl = created.upload_url || `${NOTION_API}/file_uploads/${created.id}/send`;

  // (b) 바이너리 전송 — multipart/form-data.
  //     Content-Type 은 FormData 가 boundary 와 함께 자동 설정하므로 직접 넣지 않는다.
  const form = new FormData();
  form.append('file', new Blob([buffer], { type: mime }), filename);

  const sent = await notionFetch(uploadUrl, {
      method: 'POST',
      version: FILE_UPLOAD_VERSION,
      body: form,
  });

  // ↓ 추가: 응답 상태 확인
  console.log(`    [send] status=${sent.status}, id=${sent.id}`);
  if (sent.status && sent.status !== 'uploaded') {
    throw new Error(`업로드 미완료 status=${sent.status}: ${filename}`);
  }

  const id = sent.id || created.id;
  if (!id) throw new Error(`file_upload id 를 받지 못함: ${filename}`);
  return id;
}

// ────────────────────────────────────────────────────────────────────────────
// 3) 블록 트리 후처리
// ────────────────────────────────────────────────────────────────────────────

// 이미지 블록을 file_upload 참조로 교체.
// martian 은 로컬/원격 구분 없이 image.external.url 로 변환하므로,
// url 이 우리가 업로드한 로컬 경로와 매칭되면 file_upload 로 바꾼다.
function replaceImageBlocks(blocks, uploadMap) {
  for (const block of blocks) {
    if (block.type === 'image' && block.image) {
      const url =
        (block.image.external && block.image.external.url) ||
        (block.image.file && block.image.file.url) ||
        '';
      const id = lookupUpload(url, uploadMap);
      if (id) {
        block.image = { type: 'file_upload', file_upload: { id } };
      }
    }
    // 자식 블록 재귀
    if (block[block.type] && Array.isArray(block[block.type].children)) {
      replaceImageBlocks(block[block.type].children, uploadMap);
    }
  }
  return blocks;
}

// url 을 uploadMap 키와 느슨하게 매칭 (정확/basename 양쪽 시도).
function lookupUpload(url, uploadMap) {
  if (!url) return null;
  if (uploadMap.has(url)) return uploadMap.get(url);
  const base = path.basename(url.replace(/^<|>$/g, '').trim());
  for (const [key, id] of uploadMap) {
    if (path.basename(key) === base) return id;
  }
  return null;
}

// rich_text 의 text.content 가 2000자를 넘으면 여러 element 로 분할.
function splitLongRichText(blocks) {
  for (const block of blocks) {
    const payload = block[block.type];
    if (payload && Array.isArray(payload.rich_text)) {
      payload.rich_text = payload.rich_text.flatMap(splitRichTextItem);
    }
    if (payload && Array.isArray(payload.children)) {
      splitLongRichText(payload.children);
    }
  }
  return blocks;
}

function splitRichTextItem(item) {
  const content = item && item.text && typeof item.text.content === 'string'
    ? item.text.content
    : null;
  if (content === null || content.length <= RICH_TEXT_LIMIT) return [item];

  const parts = [];
  for (let i = 0; i < content.length; i += RICH_TEXT_LIMIT) {
    const chunk = content.slice(i, i + RICH_TEXT_LIMIT);
    parts.push({
      ...item,
      text: { ...item.text, content: chunk },
    });
  }
  return parts;
}

// ────────────────────────────────────────────────────────────────────────────
// 4) 페이지 생성 + children 삽입
//    children 은 1회 100개 제한이므로 첫 100개는 페이지 생성 시,
//    나머지는 PATCH /v1/blocks/{id}/children 로 이어붙인다.
// ────────────────────────────────────────────────────────────────────────────
async function createPage(parentPageId, title, blocks) {
  const first = blocks.slice(0, CHILDREN_PER_REQUEST);
  const rest = blocks.slice(CHILDREN_PER_REQUEST);

  const page = await notionFetch(`${NOTION_API}/pages`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({
      parent: { type: 'page_id', page_id: parentPageId },
      properties: {
        title: {
          title: [{ type: 'text', text: { content: title } }],
        },
      },
      children: first,
    }),
  });

  for (let i = 0; i < rest.length; i += CHILDREN_PER_REQUEST) {
    const chunk = rest.slice(i, i + CHILDREN_PER_REQUEST);
    await notionFetch(`${NOTION_API}/blocks/${page.id}/children`, {
      method: 'PATCH',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ children: chunk }),
    });
  }

  return page;
}

// ────────────────────────────────────────────────────────────────────────────
// 페이지 제목: 첫 H1, 없으면 파일명
// ────────────────────────────────────────────────────────────────────────────
function deriveTitle(md, mdPath) {
  const h1 = md.match(/^\s*#\s+(.+?)\s*$/m);
  if (h1) return h1[1].trim();
  return path.basename(mdPath, path.extname(mdPath));
}

// ────────────────────────────────────────────────────────────────────────────
// main
// ────────────────────────────────────────────────────────────────────────────
async function main() {
  const [, , mdPathArg, parentPageId] = process.argv;

  if (!mdPathArg || !parentPageId) {
    console.error('사용법: NOTION_TOKEN=xxx node upload-to-notion.js <md경로> <부모페이지ID>');
    process.exit(1);
  }
  if (!NOTION_TOKEN) {
    console.error('환경변수 NOTION_TOKEN 이 필요합니다.');
    process.exit(1);
  }

  const mdPath = path.resolve(mdPathArg);
  if (!fs.existsSync(mdPath)) {
    console.error(`md 파일을 찾을 수 없습니다: ${mdPath}`);
    process.exit(1);
  }

  const md = fs.readFileSync(mdPath, 'utf8');
  const mdDir = path.dirname(mdPath);

  // 1) 로컬 이미지 경로 추출
  const localImages = extractLocalImagePaths(md);
  console.log(`로컬 이미지 ${localImages.length}개 발견`);

  // 2) 업로드 → (원본 md 표기 경로) → file_upload id 맵
  const uploadMap = new Map();
  for (const rel of localImages) {
    const abs = path.resolve(mdDir, rel);
    if (!fs.existsSync(abs)) {
      console.warn(`  ⚠ 이미지 없음, 건너뜀: ${rel}`);
      continue;
    }
    process.stdout.write(`  ↑ 업로드: ${rel} … `);
    const id = await uploadImage(abs);
    uploadMap.set(rel, id);
    console.log(`ok (${id})`);
  }

  // 3) md 를 이미지 기준으로 분할 변환
  //    martian 은 로컬 상대경로(및 비표준 도메인) 이미지를 드롭하므로,
  //    이미지 위치에서는 martian 을 거치지 않고 file_upload 블록을 직접 만들어 끼운다.
  //    이미지 사이의 텍스트 조각만 markdownToBlocks 로 변환한다.
  let blocks = [];

  const imgPattern =
    /!\[[^\]]*\]\(\s*<?([^)\s>]+)>?(?:\s+["'][^"']*["'])?\s*\)|<img\b[^>]*\bsrc\s*=\s*["']([^"']+)["'][^>]*>/g;

  let lastIndex = 0;
  let match;
  while ((match = imgPattern.exec(md)) !== null) {
    const imgPath = (match[1] || match[2] || '').trim();
    const before = md.slice(lastIndex, match.index);

    // 이미지 앞쪽 텍스트 → martian 블록
    if (before.trim()) {
      blocks.push(...markdownToBlocks(before));
    }

    // 이미지 → file_upload 블록 (uploadMap 에서 정확/basename 매칭)
    const id = lookupUpload(imgPath, uploadMap);
    if (id) {
      blocks.push({
        object: 'block',
        type: 'image',
        image: { type: 'file_upload', file_upload: { id } },
      });
    } else {
      console.warn(`  ⚠ 매칭 실패, 이미지 건너뜀: ${imgPath}`);
    }

    lastIndex = imgPattern.lastIndex;
  }

  // 마지막 이미지 뒤 남은 텍스트
  const tail = md.slice(lastIndex);
  if (tail.trim()) {
    blocks.push(...markdownToBlocks(tail));
  }

  console.log(
    `변환된 이미지 블록: ${blocks.filter((b) => b.type === 'image').length}개 (업로드: ${uploadMap.size}개)`
  );

  // 4) 2000자 초과 rich_text 분할
  blocks = splitLongRichText(blocks);

  // 5) 페이지 생성
  const title = deriveTitle(md, mdPath);
  console.log(`페이지 생성: "${title}" (블록 ${blocks.length}개)`);
  const page = await createPage(parentPageId, title, blocks);

  console.log('\n완료 ✅');
  console.log(`page id : ${page.id}`);
  if (page.url) console.log(`url     : ${page.url}`);
}
main().catch((err) => {
  console.error('\n실패 ❌');
  console.error(err.message);
  process.exit(1);
});
